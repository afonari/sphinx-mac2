// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_EIGENSYSTEM_H_
#define _SX_EIGENSYSTEM_H_

#include <SxMath.h>
#include <SxVector.h>
#include <SxSymMatrix.h>

/** \brief Container for eigenvalues/eigenvectors of a
           symmetric/Hermitean square matrix

    \author Christoph Freysoldt, freysoldt@mpie.de */
template<class T>
class SX_EXPORT_MATH SxSymEigensystem
{
   public:
      /// The eigenvalues
      SxVector<typename SxTypeMapper<T>::TReal> vals;
      /// The right eigenvectors
      SxVector<T> vecs;

      /// compute routine
      void compute (SxVector<T> &&mat, EIGCMD cmd, bool sorting = true);

      /// compute routine (packed matrix)
      void compute (SxVecRef<T,PackedSymMatrix> &&mat, EIGCMD cmd, bool sorting = true);

      /// Empty constructor (use compute later...)
      SxSymEigensystem () {}

      /// Constructor
      SxSymEigensystem (const SxVector<T> &mat, EIGCMD cmd = All, bool sorting = true)
      {
         compute (SxVector<T> (mat), cmd, sorting);
      }

      /// Constructor from about-to-die matrix
      SxSymEigensystem (SxVector<T> &&mat, EIGCMD cmd = All, bool sorting = true)
      {
         compute (std::move(mat), cmd, sorting);
      }

      /// Constructor: compute eigenvalues and eigenvectors
      SxSymEigensystem (const SxVector<T> &mat, bool sorting)
      {
         compute (SxVector<T> (mat), All, sorting);
      }

      /// Constructor (packed matrix)
      SxSymEigensystem (const SxSymMatrix<T> &mat, EIGCMD cmd = All, bool sorting = true)
      {
         compute (SxSymMatrix<T> (mat), cmd, sorting);
      }

      /// Constructor from about-to-die packed matrix
      SxSymEigensystem (SxSymMatrix<T> &&mat, EIGCMD cmd = All, bool sorting = true)
      {
         compute (std::move(mat), cmd, sorting);
      }

      void print (bool vectorForm=false) const  {
         SX_CHECK (vals.getSize() > 0 && vecs.getSize() > 0,
                   vals.getSize(),       vecs.getSize());
         //   SX_CHECK (vals.getSize() == vecs.getSize(),
         //               vals.getSize(),   vecs.getSize());

         printf ("Eigenvalues:\n");
         vals.print (vectorForm);
         printf ("Eigenvectors:\n");
         vecs.print (vectorForm);
      }
};

template<class T, class T2>
void sxeig_reorder (SxVector<T> *vecs, SxVector<T2> *vals, SxArray<ssize_t> &&sortIdx)
{
   SX_CHECK (vals);
   SX_CHECK (vals->getSize () == sortIdx.getSize (), vals->getSize (), sortIdx.getSize ());
   SX_CHECK (!vecs || vecs->getNCols () == sortIdx.getSize (),
             vecs->getNCols (), sortIdx.getSize ());
   SxVector<T> tmpCol;
   // --- sort in place
   for (ssize_t i = 0; i < sortIdx.getSize (); i++)  {
      if (sortIdx(i) == i) continue;
      // temp. copy value/vector at current position i
      if (vecs) tmpCol = vecs->colRef (i);
      T2 tmpVal = (*vals)(i);
      // now reorder data in a cyclic permutation
      ssize_t j = sortIdx(i), jj = i;
      while (j != i)  {
         (*vals)(jj) = (*vals)(j);
         if (vecs) vecs->colRef (jj) = vecs->colRef(j);
         sortIdx(jj) = jj;
         jj = j;
         j = sortIdx(j);
      }
      (*vals)(jj) = tmpVal;
      sortIdx(jj) = jj;
      if (vecs) vecs->colRef (jj) = tmpCol;
   }
}

template<class T>
void sxeig_GramSchmidt (SxVector<T> *vecs,
                        const SxVector<typename SxTypeMapper<T>::TReal> &vals,
                        bool sorted)
{
   SX_CHECK (vecs);
   SX_CHECK (vecs->getNCols () == vals.getSize (),
             vecs->getNCols (), vals.getSize ());
   // --- Gram-Schmidt orthogonalization of degenerate eigenvectors
   SxVecRef<T> colI;
   for (ssize_t i = 0; i < vals.getSize (); ++i)  {
      bool renormalize = false;
      for (ssize_t j = i - 1; j >= 0; --j)  {
         if (fabs(vals(i) - vals(j)) < 1e-15)  {
            if (!renormalize) colI = vecs->colRef (i);
            SxVecRef<T> colJ = vecs->colRef (j);
            colI.plus_assign_ax (-dot(colJ, colI), colJ);
            renormalize = true;
         } else if (sorted) {
            break;
         }
      }
      if (renormalize) {
         colI.normalize ();
         colI.unref ();
      }
   }
}
template<class T>
void SxSymEigensystem<T>::compute (SxVector<T> &&mat, EIGCMD cmd, bool sorting)
{
   SX_CHECK (cmd == All || cmd == ValuesOnly, cmd);
   SX_CHECK (mat.getSize () > 0);
   SX_CHECK (mat.isHermitian ());
   int n = int(mat.getNRows ());
   SX_CHECK (mat.getNRows () == n, mat.getNRows ());

   vals.resize (n);
   if (cmd == All)
      vecs.reformat (n, n);
   else
      vecs.resize (0);

   matSymEigensolver (vals.elements, vecs.elements, mat.elements, n, cmd);

   if (sorting)
      sxeig_reorder (cmd == All ? &vecs : NULL, &vals, vals.getSortIdx ());
   if (cmd == All)
      sxeig_GramSchmidt (&vecs, vals, sorting);
}

template<class T>
void SxSymEigensystem<T>::compute (SxVecRef<T,PackedSymMatrix> &&mat, EIGCMD cmd, bool sorting)
{
   SX_CHECK (cmd == All || cmd == ValuesOnly, cmd);
   SX_CHECK (mat.getSize () > 0);
   int n = int(mat.getNRows ());
   SX_CHECK (mat.getNRows () == n, mat.getNRows ());

   vals.resize (n);
   if (cmd == All)
      vecs.reformat (n, n);
   else
      vecs.resize (0);

   matEigensolverTri (vals.elements, vecs.elements, mat.elements, n, UpperRight, cmd);

   if (sorting)
      sxeig_reorder (cmd == All ? &vecs : NULL, &vals, vals.getSortIdx ());
   if (cmd == All)
      sxeig_GramSchmidt (&vecs, vals, sorting);
}

/// Get eigenvalues (symmetric matrix)
template<class T>
SxVector<typename SxTypeMapper<T>::TReal> symEigenvalues (const SxVector<T> &mat)  {
   return std::move(SxSymEigensystem<T> (mat, ValuesOnly).vals);
}

/// Get eigenvalues (about-to-die symmetric matrix)
template<class T>
SxVector<T> symEigenvalues (SxVector<T> &&mat)  {
   return std::move(SxSymEigensystem<T> (std::move(mat), ValuesOnly).vals);
}

/// Get eigenvalues (symmetric packed matrix)
template<class T>
SxVector<typename SxTypeMapper<T>::TReal> symEigenvalues (const SxSymMatrix<T> &mat)  {
   return std::move(SxSymEigensystem<T> (mat, ValuesOnly).vals);
}

/// Get eigenvalues (about-to-die symmetric packed matrix)
template<class T>
SxVector<T> symEigenvalues (SxSymMatrix<T> &&mat)  {
   return std::move(SxSymEigensystem<T> (std::move(mat), ValuesOnly).vals);
}

/** \brief Container for eigenvalues/eigenvectors of a general square matrix

    \author Christoph Freysoldt, freysoldt@mpie.de */
template<class T>
class SX_EXPORT_MATH SxEigensystem
{
   public:
      /// The eigenvalues
      SxVector<typename SxTypeMapper<T>::TCmplx> vals;
      /** \brief The right eigenvectors
          \note For real matrices and complex eigenvalue pairs (a +- bi), the
                complex eigenvector pair is stored with separate real/imag
                values (cf. LAPACK docu)
        */
      SxVector<T> vecs;

      /// compute routine
      void compute (SxVector<T> &&mat, EIGCMD cmd, bool sorting);

      /// Constructor
      SxEigensystem (const SxVector<T> &mat, EIGCMD cmd = All, bool sorting = false)
      {
         compute (SxVector<T> (mat), cmd, sorting);
      }

      /// Constructor from about-to-die matrix
      SxEigensystem (SxVector<T> &&mat, EIGCMD cmd = All, bool sorting = false)
      {
         compute (std::move(mat), cmd, sorting);
      }

      void print (bool vectorForm=false) const  {
         SX_CHECK (vals.getSize() > 0 && vecs.getSize() > 0,
                   vals.getSize(),       vecs.getSize());
         //   SX_CHECK (vals.getSize() == vecs.getSize(),
         //               vals.getSize(),   vecs.getSize());

         printf ("Eigenvalues:\n");
         vals.print (vectorForm);
         printf ("Eigenvectors:\n");
         vecs.print (vectorForm);
      }
};

template<class T>
void SxEigensystem<T>::compute (SxVector<T> &&mat, EIGCMD cmd, bool sorting)
{
   SX_CHECK (cmd == All || cmd == ValuesOnly, cmd);
   SX_CHECK (mat.getSize () > 0);
   int n = int(mat.getNRows ());
   SX_CHECK (mat.getNRows () == n, mat.getNRows ());

   vals.resize (n);
   if (cmd == All)
      vecs.reformat (n, n);
   else
      vecs.resize (0);

   matEigensolver (vals.elements, vecs.elements, mat.elements, n, cmd);
   if (sorting)  {
      SxVector<typename SxTypeMapper<T>::TReal> valsRe = vals.real ();
      sxeig_reorder (cmd == All ? &vecs : NULL, &vals, valsRe.getSortIdx ());
   }
}

template<class T>
SxVector<T> eigenvalues (const SxVector<T> &mat, bool sorting = true)  {
   return std::move(SxEigensystem<T> (mat, ValuesOnly, sorting).vals);
}

template<class T>
SxVector<T> eigenvalues (SxVector<T> &&mat, bool sorting = true)  {
   return std::move(SxEigensystem<T> (std::move(mat), ValuesOnly, sorting).vals);
}

#endif /* _SX_EIGENSYSTEM_H_ */
