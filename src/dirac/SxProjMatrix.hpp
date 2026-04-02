// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institut fuer Eisenforschung GmbH
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxrepo.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

template<class T>
SxVector<T>
SxProjMatrix<T>::getProjectedSingle (int iStart, int iEnd,
                                     const SxVecRef<T> &projBlock,
                                     const SxVecRef<T> &right) const
{
   SX_CHECK (iEnd - iStart + 1 == projBlock.getNCols (),
             iEnd, iStart, projBlock.getNCols ());
   // return (projBlock.adjoint () ^ right); // XPress
   return (right.conj () ^ projBlock).conj ()
          .reshape (projBlock.getNCols ());
}

template<class T>
SxVector<T>
SxProjMatrix<T>::SaveProjections::getProjectedSingle (int iStart, int iEnd,
                                     const SxVecRef<T> &projBlock,
                                     const SxVecRef<T> &right) const
{
   SX_CHECK (iStart + projBlock.getNCols () == iEnd + 1,
             iStart, projBlock.getNCols (), iEnd);
   if (cacheUsed)  {
      return savedProjections.getRef(projBlock.getNCols (), iStart);
   }
   SxVector<T> res = SxProjMatrix<T>::getProjectedSingle (iStart, iEnd,
                                                          projBlock, right);
   if (useCache)
      savedProjections.getRef(projBlock.getNCols (), iStart) = res;
   return res;
}

template<class T>
T
SxProjMatrix<T>::getProjectedSingle (int,
                                     const SxVecRef<T> &projector,
                                     const SxVecRef<T> &right) const
{
   return dot (projector, right);
}


template<class T>
T
SxProjMatrix<T>::SaveProjections::getProjectedSingle (
      int iProj,
      const SxVecRef<T> &projector,
      const SxVecRef<T> &right) const
{
   if (cacheUsed)  {
      return savedProjections(iProj);
   } else {
      // perform projection
      T proj = dot (projector, right);
      if (useCache) savedProjections(iProj) = proj;
      return proj;
   }
}
      
template <class T>
int
SxProjMatrix<T>::SaveProjections::prepare (const SxVecRef<T> &right) const
{
   int nR = (int)right.getNCols ();
   if (nR == 0 && right.getSize () == nElements) nR = 1;
   if (useCache)  {
      savedProjections = SxVector<T> (nProj * nR);
      savedProjections.reshape (nProj, nR);
      cacheUsed = false;
   } else if (cacheUsed) {
      nR = (int)savedProjections.getNCols ();
   }
   return nR;
}

template <class T>
void SxProjMatrix<T>::SaveProjections::finalize () const
{
   if (useCache) cacheUsed = true;
}

template<class T>
void SxProjMatrix<T>::applyAndAddSingle   (      SxVecRef<T> *left,
                                           const SxVecRef<T> &right) const
{
   SX_CHECK (left);
   SX_CHECK (left->getSize () == nElements, left->getSize (), nElements);
   
   // hook in for derived classes
   prepare (right);

   int ip, iBlock = 0;
   SxVector<T> projBlock;
   projBlock.reformat (nElements, blocksize);

   for (ip=0; ip < nProj; ++ip)  {
      
      // collect projectors
      SxVecRef<T> projCol = projBlock.colRef(iBlock++);
      getProjector(ip, &projCol);

      // --- when block is full
      if (iBlock == projBlock.getNCols () /* == current blocksize */)  {
         if (iBlock > 1)  {
            SxVector<T> projected
               = getProjectedSingle (ip-iBlock+1, ip, projBlock, right);
            
            if (useFactors ())  {
               // multiply with inner factor
               typename SxVector<T>::Iterator pIt = projected.begin ();
               for (int i = ip - iBlock + 1; i <= ip; ++i)
                  *pIt++ *= getFactor (i);
            }

            // add to left side
            *left += projBlock ^ projected;
         } else {
            T proj = getProjectedSingle (ip, projBlock, right);
            // multiply with inner factor
            proj *= getFactor (ip);

            // add to left side
            left->plus_assign_ax (proj, projBlock);
         }
         
         iBlock = 0;
         // reformat next block
         if (blocksize + ip >= nProj && ip + 1 < nProj)
            projBlock.reformat (nElements, nProj - ip - 1);
      }
   }
   // hook in for derived classes
   finalize ();
}

template<class T>
SxVector<T>
SxProjMatrix<T>::getProjectedMultiple (int iStart, int iEnd,
                                     const SxVecRef<T> &projBlock,
                                     const SxVecRef<T> &right) const
{
   SX_CHECK (iEnd - iStart + 1 == projBlock.getNCols (),
             iEnd, iStart, projBlock.getNCols ());
   //return (projBlock.adjoint () ^ right);
   // note: "right" may have more rows than we use
   return projBlock.overlap (right, projBlock.getNRows ());
}

template<class T>
SxVector<T>
SxProjMatrix<T>::SaveProjections::getProjectedMultiple (int iStart, int iEnd,
                                     const SxVecRef<T> &projBlock,
                                     const SxVecRef<T> &right) const
{
   SxVector<T> res;
   if (cacheUsed)  {
      int nCols = (int)savedProjections.getNCols ();
      int nRows = iEnd - iStart + 1;
      res.reformat (nRows, nCols);
      for (int ic = 0; ic < nCols; ++ic)  {
         int offset = nProj * ic + iStart;
         res.colRef(ic) <<= savedProjections.getRef (nRows, offset);
      }
   } else {
      res = SxProjMatrix<T>::getProjectedMultiple (iStart, iEnd, projBlock,
                                                 right);
      if (useCache)  {
         int nCols = (int)res.getNCols ();
         int nRows = iEnd - iStart + 1;
         for (int ic = 0; ic < nCols; ++ic)  {
            int offset = nProj * ic + iStart;
            savedProjections.getRef (nRows, offset) = res.colRef(ic);
         }
      }
   }
   return res;
}

template<class T>
SxVector<T>
SxProjMatrix<T>::getProjectedMultiple (int,
                                     const SxVecRef<T> &projector,
                                     const SxVecRef<T> &right) const
{
   //return projector.conj ().reshape (1, nElements) ^ right;
   // note: "right" may have more rows than we use
   SxVector<T> res = right.overlap (projector, projector.getNRows ());
   // form conjugate in place
   for (int i = 0; i < res.getSize (); i++)
      res(i).im = -res(i).im;
   return res;
}


template<class T>
SxVector<T>
SxProjMatrix<T>::SaveProjections::getProjectedMultiple (
      int iProj,
      const SxVecRef<T> &projector,
      const SxVecRef<T> &right) const
{
   SxVector<T> res;
   if (cacheUsed)  {
      // get a single row
      res = savedProjections.rowRef (iProj);
   } else {
      // perform projection
      res = SxProjMatrix<T>::getProjectedMultiple(iProj, projector, right);
      if (useCache) savedProjections.rowRef (iProj) = res;
   }
   return res;
}
      
template<class T>
void SxProjMatrix<T>::applyAndAddMultiple (      SxVecRef<T> *left,
                                           const SxVecRef<T> &right) const
{
   SX_CHECK (left);
   
   int nCols = prepare (right);

   SX_CHECK (nCols > 1);
   SX_CHECK (left->getNRows () == nElements, left->getNRows (), nElements);
   SX_CHECK (left->getNCols () == nCols, left->getNCols (), nCols);

   int ip, iBlock = 0;
   SxVector<T> projBlock (nElements, blocksize);

   for (ip=0; ip < nProj; ++ip)  {
      
      // collect projectors
      SxVecRef<T> projCol = projBlock.colRef(iBlock++);
      getProjector(ip, &projCol);

      // --- when block is full
      if (iBlock == projBlock.getNCols () /* == current blocksize */)  {
         if (iBlock > 1)  {
            SxVector<T> projected
               = getProjectedMultiple (ip - iBlock + 1, ip, projBlock, right);
            
            // get inner factor
            if (useFactors ())  {
               SxVector<T> diag(iBlock);
               typename SxVector<T>::Iterator diagIt = diag.begin ();
               for (int i = ip - iBlock + 1; i <= ip; ++i)
                  *diagIt++ = getFactor (i);
               // multiply with inner factor
               for (int ic = 0; ic < nCols; ++ic)
                  projected.colRef(ic) *= diag;
            }

            // add to left side
            *left += projBlock ^ projected;
         } else {
            SxVector<T> projected = getProjectedMultiple (ip, projBlock, right);
            
            // multiply with inner factor
            T factor = getFactor (ip);

            // add to left side
            typename SxVector<T>::Iterator pIt = projected.begin ();
            for (int ic = 0; ic < nCols; ++ic, ++pIt)
               left->colRef(ic).plus_assign_ax(factor * *pIt,
                                               projBlock);
         }
         
         iBlock = 0;
         // reformat next block
         if (blocksize + ip >= nProj && ip + 1 < nProj)
            projBlock.reformat (nElements, nProj - ip - 1);
      }
   }
   // hook in for derived classes
   finalize ();
}

template<class T>
void SxProjMatrix<T>::projectionKernel (SxVector<T> *res,
                                        const SxVecRef<T> &left,
                                        const SxVecRef<T> &right) const
{
   int nL = (int)left.getNCols ();
   if (nL == 0 && left.getSize () == nElements) nL = 1;
   SX_CHECK (nL == 0 || res, nL);
   int nR = prepare (right);

   if (res)  {
      if (res->getSize () > 0 && res->getNRows () == nL && res->getNCols () == nR)  {
         res->set (0.);
      } else {
         *res = SxVector<T> ();
      }
   }

   int ip, iBlock = 0;
   SxVector<T> projBlock (nElements, blocksize);

   for (ip=0; ip < nProj; ++ip)  {
      
      // collect projectors
      SxVecRef<T> projCol = projBlock.colRef(iBlock++);
      getProjector(ip, &projCol);

      // --- when block is full
      if (iBlock == projBlock.getNCols () /* == current blocksize */)  {
         if (iBlock > 1)  {
            // --- perform projection
            SxVector<T> projR = getProjectedMultiple (ip-iBlock+1,ip,projBlock,right);
            projR.reshape (iBlock, nR);

            if (nL > 0)  {
               SxVecRef<T> projL;
               
               if (&left == &right)
                  projL = projR;
               else
                  //projL = projBlock.adjoint () ^ left;
                  projL = projBlock.overlap (left);
               // may not be the case for nL == 1
               projL.reshape (iBlock, nL);
               
               // --- calculate contribution to matrix

               if (useFactors ())  {
                  if (&left == &right) SX_EXIT; // CF 2020-05-23: projL references projR
                  // get inner factor
                  SxVector<T> diag(iBlock);
                  typename SxVector<T>::Iterator diagIt = diag.begin ();
                  for (int i = ip - iBlock + 1; i <= ip; ++i)
                     *diagIt++ = getFactor (i);
                  // multiply with inner factor
                  for (int ic = 0; ic < nR; ++ic)
                     projR.colRef(ic) *= diag;
               }

               // get contribution of this block
               SxVector<T> LR = projL.overlap (projR);

               // add to result
               if (res->getSize () == 0)
                  *res = LR;
               else
                  *res += LR;
            }
         } else {
            // --- perform projection
            SxVector<T> projR = getProjectedMultiple (ip, projBlock, right);
            projR.reshape (1, nR);
            if (nL > 0)  {
               SxVecRef<T> projL;
               if (&left == &right)
                  projL = projR.adjoint ();
               else
                  //projL = (projBlock.adjoint () ^ left).adjoint ();
                  projL = (projBlock.overlap (left)).adjoint ();
               projL.reshape (nL, 1);
               
               // --- calculate contribution to matrix
               
               // get inner factor
               T factor = getFactor (ip);

               if (res->getSize () == 0)  {
                  res->reformat (nL, nR);
                  res->set (0.);
               }

               // add contribution of this block to result
               if (nL == 1)  {
                  res->plus_assign_ax(projL(0) * factor, projR);
               } else if (nR == 1)  {
                  res->plus_assign_ax(projR(0) * factor, projL);
               } else  {
                  // result += factor * (projL ^ projR);
                  res->plus_assign_ax (factor, projL ^ projR);
               }
            }
         }
         
         iBlock = 0;
         // reformat next block
         if (blocksize + ip >= nProj && ip + 1 < nProj)
            projBlock.reformat (nElements, nProj - ip - 1);
      }
   }
   finalize ();
}


template <class T>
void SxProjMatrix<T>::getProjector(int i, SxVecRef<T> *target) const
{
   SX_CHECK (target);
   (*target) <<= getProjector (i);
}


template <class T>
void SxProjMatrix<T>::SaveProjections::clearCache ()
{ 
   savedProjections = SxVector<T> ();
   cacheUsed = false;
}


template <class T>
SxVector<T>
SxProjMatrix<T>::SaveProjections::getProjection(const SxVecRef<T> &right)
{
   SX_CHECK (right.getSize () > 0);
   SX_CHECK(right.getNRows () == nElements, right.getNRows (), nElements);
   useCache = true;
   clearCache ();
   projectionKernel (NULL, SxVector<T> (), right);
   SX_CHECK (cacheUsed);
   return savedProjections;
}

template <class T>
SxVector<T>
SxProjMatrix<T>::SaveProjections::getProjectionFromExtended(const SxVecRef<T> &right)
{
   SX_CHECK (right.getSize () > 0);
   SX_CHECK(right.getNRows () >= nElements, right.getNRows (), nElements);
   useCache = true;
   clearCache ();
   projectionKernel (NULL, SxVector<T> (), right);
   SX_CHECK (cacheUsed);
   return savedProjections;
}

template<class T>
SxVector<T>
SxProjMatrix<T>::SaveProjections::gradient (const SxVecRef<T> &proj)
{
   SX_CHECK (proj.getSize () > 0 || cacheUsed);
   if (proj.getSize () > 0)  {
      SX_CHECK (proj.getNRows () == nProj, proj.getNRows (), nProj);
      savedProjections = proj;
      cacheUsed = true;
   }
   useCache = false;
   int nCols = (int)savedProjections.getNCols ();
   if (nCols == 0) nCols = 1;
   SxVector<T> res (nElements, nCols);
   res.set (0.);

   if (nCols == 1)
      applyAndAddSingle (&res, SxVector<PrecCoeffG> ());
   else
      applyAndAddMultiple (&res, SxVector<PrecCoeffG> ());

   return res;
}
