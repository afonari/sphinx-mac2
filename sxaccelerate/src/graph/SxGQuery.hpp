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


// --- returns the matched selections of nodes' Idx.
template<class N,class E,template<class,bool> class GS>
SxGQuery::SelSet
SxGQuery::matchAll (const SxPtr<SxGraph<N,E,GS> > &g,
                    const typename SxGraph<N,E,GS>::Iterator &gIt)
{
   SX_TRACE ();
   selections = SxGQuery::SelSet::create ();
   SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();
   bool res = false, finalRes = false;
   typename SxGraph<N,E,GS>::ConstIterator tIt;
   // --- list for visited pattern nodes overall
   SxGQuery::Selection gVisited = SxGQuery::Selection::create ();
   SxGQuery::Selection visited = SxGQuery::Selection::create ();
   SxGQuery::Selection fullVisited = SxGQuery::Selection::create ();

   auto it = gIt;

   // --- iterate over all pattern nodes irrespective
   //     of edges to cover all connected components.
   auto patternIt = patternGraph.begin (sx::Undefined, (*patternGraph.begin ()));
   while (patternIt.isValid ())  {
      while (it.isValid ())  {
         SX_CHECK (selections.getPtr ());
         SX_CHECK (tmpSels.getPtr ());
         visited = SxGQuery::Selection::create ();
         // --- start a new iterator in order to ignore the
         //     already visited nodes of 'it'.
         tIt = g->begin (*it);
         res = matchAllRecursive<N,E,GS> (patternGraph.begin (*patternIt),
                                          tIt, tmpSels, visited);
         if (res == true)  {
            // --- at least one match found during this loop
            finalRes = true;
            // --- store the full matched set of pattern nodes
            fullVisited = visited;
            selections->append (std::move(*tmpSels));
         }
         ++it;
      }
      if (!finalRes)  {
         return SxGQuery::SelSet::create ();
      }
      gVisited->append (*fullVisited);
      finalRes = false;

      it = gIt;

      // --- skip the visited pattern nodes
      while (patternIt.isValid () && gVisited->contains (patternIt->getId ())) {
         ++patternIt;
      }
   }
   return selections;
}

template<class N,class E,template<class,bool> class GS>
SxGQuery::SelSet SxGQuery::matchAll (const SxPtr<SxGraph<N,E,GS> > &g)
{
   SX_TRACE ();
   return matchAll<N,E,GS> (g, g->begin (sx::Undefined, (*g->begin ())));
}

template<class N,class E,template<class,bool> class GS>
SxGQuery::Selection
SxGQuery::match (const SxPtr<SxGraph<N,E,GS> > &g,
                 const typename SxGraph<N,E,GS>::Iterator &gIt)
{
   SX_TRACE ();
   bool res = false;
   SxGQuery::Selection sel = SxGQuery::Selection::create ();
   SxGQuery::Selection visited = SxGQuery::Selection::create ();
   auto it = gIt;
   typename SxGraph<N,E,GS>::ConstIterator tIt;

   auto patternIt = patternGraph.begin (sx::Undefined, (*patternGraph.begin ()));
   while (patternIt.isValid ())  {
      while (it.isValid ())  {
         tIt = g->begin (*it);
         res = matchRecursive<N,E,GS> (patternGraph.begin (*patternIt),
                                       tIt, sel, visited);
         if (res == true)  break;
         ++it;
      }
      if (!res)  {
         sel->removeAll ();
         break;
      }
      res = false;
      it = gIt;
      // --- skip the visited pattern nodes
      while (patternIt.isValid () && visited->contains (patternIt->getId ()))  {
         ++patternIt;
      }
   }

   return sel;
}

template<class N,class E,template<class,bool> class GS>
SxGQuery::Selection SxGQuery::match (const SxPtr<SxGraph<N,E,GS> > &g)
{
   SX_TRACE ();
   return match<N,E,GS> (g, g->begin (sx::Undefined, (*g->begin ())));
}

// ---------------------------------------------------------------------------
// protected functions

void SxGQuery::resizeSelections (const SxGQuery::SelSet &sels,
                                 ssize_t setSize,
                                 ssize_t selSize) const
{
   SX_TRACE ();
   SX_CHECK (sels.getPtr ());

   sels->resize (setSize);

   for (auto it = sels->begin (); it != sels->end (); ++it)  {
      (*it)->resize (selSize);
   }
}

SxGQuery::SelSet SxGQuery::crossSelections (const SxGQuery::SelSet &sels1,
                                            const SxGQuery::SelSet &sels2) const
{
   SX_TRACE ();
   SX_CHECK (sels1.getPtr ());
   SX_CHECK (sels2.getPtr ());
   SxGQuery::SelSet res = SxGQuery::SelSet::create ();

   if (sels1->getSize () == 0)  {
      for (auto it2 = sels2->begin (); it2 != sels2->end (); ++it2)  {
         SxGQuery::Selection ptr = SxGQuery::Selection::create ();
         ptr->append ((*(*it2)));
         res->append (ptr);
      }
   }  else if (sels2->getSize () == 0)  {
      for (auto it1 = sels1->begin (); it1 != sels1->end (); ++it1)  {
         SxGQuery::Selection ptr = SxGQuery::Selection::create ();
         ptr->append ((*(*it1)));
         res->append (ptr);
      }
   }  else  {
      for (auto it1 = sels1->begin (); it1 != sels1->end (); ++it1)  {
         for (auto it2 = sels2->begin (); it2 != sels2->end (); ++it2)  {
            SxGQuery::Selection ptr = SxGQuery::Selection::create ();
            ptr->append ((*(*it1)));
            ptr->append ((*(*it2)));
            res->append (ptr);
         }
      }
   }

   SX_CHECK (res.getPtr ());
   return res;
}

// --- match all functions

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllInNode (const SxGraph<SxGQPattern>::ConstIterator &it,
                          const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                          const SxGQuery::SelSet &sels,
                          const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
   bool finalRes = false;
   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);
   SxGQuery::Selection currVisited = SxGQuery::Selection::create (*visited);
   for (ssize_t i = 0; i < gIt.getSizeIn (); ++i)  {
      SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();
      if(matchAllRecursive<N,E,GS> (it, gIt.in (i), tmpSels, tmpVisited))  {
         if (!finalRes)  *currVisited = *tmpVisited;
         // --- atleast one match found
         finalRes = true;
         resSet->append (std::move(*crossSelections (tmpSels, sels)) );
      }

      tmpVisited->resize (visited->getSize ());
      SX_CHECK (*tmpVisited == *visited);
   }

   if (finalRes)  {
      *sels = *resSet;
      // --- update visited list
      *visited = std::move(*currVisited);
   }
   return finalRes;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllInNodes (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                           const SxGQuery::SelSet &sels,
                           const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
   resSet= crossSelections (resSet, sels);

   SxGQuery::Selection currVisited = SxGQuery::Selection::create (*visited);

   bool res = true, hasVisited = false;

   for (ssize_t i = 0; i < it.getSizeIn (); ++i)  {
      auto inIt = it.in (i);
      if (!visited->contains (inIt->getId ()))  {
         hasVisited = true;
         res = matchAllInNode<N,E,GS> (inIt, gIt, resSet, currVisited);
         if (!res)  break;
      }
   }

   if (hasVisited && res)  {
      *sels = *resSet;
      *visited = std::move (*currVisited);
   }
   return res;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllRecursive (const SxGraph<SxGQPattern>::ConstIterator &it,
                             const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                             const SxGQuery::SelSet &sels,
                             const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   // --- if this node has already been visited
   //     still evaluate/match it so that calling function
   //     can know success/failure. But don't append
   //     it to the result.
   if (visited->contains (it->getId ()))  {
      SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();
      return it->matchAll<N,E,GS> (gIt, tmpSels);
   }

   bool res = false;
   ssize_t setSize = sels->getSize ();
   ssize_t siz = (setSize == 0)? 0 : sels->first()->getSize ();

   res = it->matchAll<N,E,GS> (gIt, sels);

   if (res)  {
      visited->append (it->getId ());
      res = matchAllInNodes<N,E,GS> (it, gIt, sels, visited);
   }

   if (res && it.getSizeOut () > 0)  {
      switch (it->getRelType ())  {
         case SxGQPattern::RelationType::OrderedDirect:
            res = matchAllOutNodesOD<N,E,GS> (it, gIt, sels, visited);
            break;
         case SxGQPattern::RelationType::OrderedIndirect:
            res = matchAllOutNodesOI<N,E,GS> (it, gIt, sels, visited);
            break;
         case SxGQPattern::RelationType::UnorderedDirect:
            res = matchAllOutNodesUD<N,E,GS> (it, gIt, sels, visited);
            break;
         case SxGQPattern::RelationType::UnorderedIndirect:
            res = matchAllOutNodesUI<N,E,GS> (it, gIt, sels, visited);
            break;
      }
   }

   if (!res)  resizeSelections (sels, setSize, siz);

   return res;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllOutNodesODIdx (const SxGraph<SxGQPattern>::ConstIterator &it,
                                 const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                 ssize_t childIdx,
                                 const SxGQuery::SelSet &sels,
                                 const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SX_CHECK (it.getSizeOut () > 0, it.getSizeOut ());

   bool res = true;
   ssize_t setSize = sels->getSize ();
   ssize_t siz = (setSize == 0)? 0 : sels->first()->getSize ();
   ssize_t outSize = gIt.getSizeOut ();

   auto currIt = gIt.out (childIdx);

   if ((outSize-childIdx) < it.getSizeOut ())  {
      return false;
   }

   for (ssize_t j = 0; j < it.getSizeOut (); ++j)  {

      res = matchAllRecursive<N,E,GS> (it.out(j), gIt.out (childIdx),
                                       sels, visited);
      if (!res)  break;
      SX_CHECK (childIdx < outSize, childIdx, outSize);
      ++childIdx;
   }

   if (!res)  resizeSelections (sels, setSize, siz);
   return res;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllOutNodesOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                              const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                              const SxGQuery::SelSet &sels,
                              const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
   bool finalRes = false;
   ssize_t outSize = gIt.getSizeOut ();
   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);
   SxGQuery::Selection currVisited = SxGQuery::Selection::create (*visited);
   for (ssize_t c = 0; c < outSize; ++c)  {
      SxGQuery::SelSet childSels = SxGQuery::SelSet::create ();

      if (matchAllOutNodesODIdx<N,E,GS> (it, gIt, c, childSels, tmpVisited))  {
         // --- store the visited list of successful match
         if (!finalRes)  *currVisited = *tmpVisited;
         finalRes = true;
         resSet->append (std::move(*crossSelections(sels, childSels)));
      }
      tmpVisited->resize (visited->getSize ());
      SX_CHECK (*tmpVisited == *visited);
   }

   if (finalRes)  {
      *sels = *resSet;
      *visited = std::move(*currVisited);
   }
   return finalRes;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllOutNodesUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                              const typename SxGraph<N,E,GS>::ConstIterator &pIt,
                              const SxGQuery::SelSet &sels,
                              const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SxGQuery::Selection prevMatches = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > prevSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   ssize_t outSize = pIt.getSizeOut ();

   SxGQuery::Selection currVisited = SxGQuery::Selection::create (*visited);
   SxGQuery::Selection tmpVisited  = SxGQuery::Selection::create (*currVisited);

   // -- list of visited pattern nodes for a particular pattern
   SxGQuery::Selection successVisited = SxGQuery::Selection::create (*visited);

   bool finalRes = false;
   bool res = false;

   SxGQuery::Selection currMatches   = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > currSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   for (ssize_t exprIdx = 0; exprIdx < it.getSizeOut (); ++exprIdx)  {

      auto chIt = it.out (exprIdx);
      finalRes = false;
      res      = false;

      *tmpVisited = *currVisited;

      for (ssize_t cIdx = 0; cIdx < outSize; ++cIdx)  {

         SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();

         res = matchAllRecursive<N,E,GS> (chIt, pIt.out (cIdx),
                                          resSet, tmpVisited);
         if (res == true)  {
            // --- store the visited list of successful match
            if (!finalRes)  *successVisited = *tmpVisited;

            // --- atleast one match found for curr expr
            finalRes = true;
            currMatches->append (cIdx);
            SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();

            if (prevMatches->getSize () > 0)  {

               auto pMIt = prevMatches->begin ();
               auto pSIt = prevSels->begin ();
               for (;pMIt != prevMatches->end (); (++pMIt,++pSIt))  {
                  tmpSels->append (std::move( *(crossSelections (*pSIt, resSet)) ) );
               }

            }  else  {
               // --- prevMatches is empty i-e. we are in first expression
               *tmpSels = *resSet;
            }

            SX_CHECK (tmpSels->getSize () > 0, tmpSels->getSize ());
            currSels->append (tmpSels);

         }

         tmpVisited->resize (currVisited->getSize ());
         SX_CHECK (*tmpVisited == *currVisited);
      }

      // --- if no match found then return
      if (finalRes == false)  {
         return false;
      }

      // --- update the visited pattern nodes uptill now
      *currVisited = std::move(*successVisited);


      *prevMatches = std::move(*currMatches);
      *prevSels = std::move(*currSels);

   }

   SxGQuery::SelSet rSet = SxGQuery::SelSet::create ();
   for (auto sIt = prevSels->begin (); sIt != prevSels->end (); ++sIt)  {
      rSet->append ( std::move( *(crossSelections (sels, *sIt)) ) );
   }

   *visited = std::move(*currVisited);

   *sels = *rSet;

   return true;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllOutNodesOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                              const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                              const SxGQuery::SelSet &sels,
                              const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SxGQuery::Selection prevMatches = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > prevSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   ssize_t cIdx = 0;
   ssize_t outSize = gIt.getSizeOut ();
   ssize_t exprsLeft = it.getSizeOut ();


   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);
   SxGQuery::Selection currVisited = SxGQuery::Selection::create (*visited);
   SxGQuery::Selection successVisited = SxGQuery::Selection::create (*visited);


   bool finalRes = false;
   bool res = false;

   SxGQuery::Selection currMatches   = SxGQuery::Selection::create ();
   SxPtr<SxList<SxGQuery::SelSet> > currSels = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   for (ssize_t exprIdx = 0; exprIdx < it.getSizeOut (); ++exprIdx)  {
      auto chIt = it.out (exprIdx);

      finalRes = false;
      res      = false;

      *tmpVisited = *currVisited;

      for (; cIdx < outSize; ++cIdx)  {
         // --- if number of siblings left to process
         //     are less than the number of expresssions
         //     left then stop
         if ((outSize-cIdx) < exprsLeft)
            break;

         SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
         res = matchAllRecursive<N,E,GS> (chIt, gIt.out (cIdx), resSet, tmpVisited);
         if (res == true)  {
            if (!finalRes)  *successVisited = *tmpVisited;
            // --- atleast one match found for curr expr
            finalRes = true;
            currMatches->append (cIdx);
            SxGQuery::SelSet tmpSels = SxGQuery::SelSet::create ();

            if (prevMatches->getSize () > 0)  {

               auto pMIt = prevMatches->begin ();
               auto pSIt = prevSels->begin ();
               for (;pMIt != prevMatches->end (); (++pMIt,++pSIt))  {

                  if (*pMIt < cIdx)  {
                     tmpSels->append (std::move( *(crossSelections (*pSIt, resSet)) ) );
                  }  else  {
                     // --- childIdx's in prevMatches are in sorted order
                     //     so once (*pMIt >= cIdx), the rest of them will
                     //     also be greater than cIdx. Due to "Ordered"
                     //     requirement those won't match. so break loop.
                     break;
                  }

               }
            }  else  {
               // --- prevMatches is empty i-e. in first expression
               *tmpSels = *resSet;
            }

            SX_CHECK (tmpSels->getSize () > 0, tmpSels->getSize ());
            currSels->append (tmpSels);
         }

         tmpVisited->resize (currVisited->getSize ());
         SX_CHECK (*tmpVisited == *currVisited);
      }

      // --- if no match found then return
      if (finalRes == false)  return false;

      *currVisited = std::move(*successVisited);

      --exprsLeft;
      *prevMatches = std::move (*currMatches);
      *prevSels = std::move (*currSels);

      // --- start next loop from firstChildIdx+1
      //     that matched in prev iteration.
      SX_CHECK (prevMatches->getSize() > 0, prevMatches->getSize ());
      cIdx = prevMatches->first ()+1;

   }

   SxGQuery::SelSet rSet = SxGQuery::SelSet::create ();
   for (auto sIt = prevSels->begin (); sIt != prevSels->end (); ++sIt)  {
      rSet->append ( std::move( *(crossSelections (sels, *sIt)) ) );
   }

   *sels = *rSet;

   *visited = std::move (*currVisited);

   return true;

}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchAllOutNodesUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                              const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                              const SxGQuery::SelSet &sels,
                              const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SX_CHECK (it.getSizeOut () > 0, it.getSizeOut ());

   SxGQuery::Selection prevMatches = SxGQuery::Selection::create ();

   SxGQuery::Selection currVisited = SxGQuery::Selection::create (*visited);
   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);
   SxGQuery::Selection successVisited = SxGQuery::Selection::create ();

   SxPtr<SxList<SxGQuery::SelSet> > prevSels
      = SxPtr<SxList<SxGQuery::SelSet> >::create ();
   ssize_t outSize = gIt.getSizeOut ();
   ssize_t cIdx    = 0;
   ssize_t exprIdx = 0;

   bool finalRes = false;
   bool res = false;
   bool res1 = false, res2 = false;

   auto exprIt = it.out (exprIdx);

   // --- go over the childList to find all matches for first expr
   for (; cIdx < outSize; ++cIdx)  {

      SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
      res = matchAllRecursive<N,E,GS> (exprIt, gIt.out (cIdx), resSet, tmpVisited);
      if (res == true)  {
         if (!finalRes)  *successVisited = *tmpVisited;
         // --- at least one match found for curr expr
         finalRes = true;
         prevMatches->append (cIdx);
         prevSels->append (resSet);
      }

      // resize back
      tmpVisited->resize (currVisited->getSize ());
      SX_CHECK (*tmpVisited == *currVisited);
   }



   if (!finalRes)  return false;

   *currVisited = std::move(*successVisited);


   exprIdx++;

   for (; exprIdx < it.getSizeOut (); ++exprIdx)  {

      exprIt = it.out (exprIdx);

      SxGQuery::Selection currMatches = SxGQuery::Selection::create ();
      SxPtr<SxList<SxGQuery::SelSet> > currSels =
                                    SxPtr<SxList<SxGQuery::SelSet> >::create ();
      finalRes = false;
      res      = false;
      auto pMIt = prevMatches->begin ();
      auto pSIt = prevSels->begin ();

      *tmpVisited = *currVisited;

      for (; pMIt != prevMatches->end (); ++pMIt, ++pSIt)  {

         cIdx = *pMIt;
         res1 = false; res2 = false;
         if ((cIdx-1) >= 0)  {
            auto prevIt = gIt.out ((cIdx-1));
            SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
            res1 = matchAllRecursive<N,E,GS> (exprIt, prevIt, resSet, tmpVisited);
            if (res1)  {
               *successVisited = *tmpVisited;
               currMatches->append ((cIdx-1));
               currSels->append ( crossSelections (*pSIt, resSet) );
            }  else  {
               tmpVisited->resize (currVisited->getSize ());
            }
            tmpVisited->resize (currVisited->getSize ());
            SX_CHECK (*tmpVisited == *currVisited);
         }

         if ((cIdx+1) < outSize)  {
            auto nextIt = gIt.out ((cIdx+1));
            SxGQuery::SelSet resSet = SxGQuery::SelSet::create ();
            res2 = matchAllRecursive<N,E,GS> (exprIt, nextIt, resSet, tmpVisited);
            if (res2)  {
               *successVisited = *tmpVisited;
               currMatches->append ((cIdx+1));
               currSels->append ( crossSelections (*pSIt, resSet) );
            }
            tmpVisited->resize (currVisited->getSize ());
            SX_CHECK (*tmpVisited == *currVisited);
         }

         res = res1 || res2;

         // --- if any one of the prevMatched child results in
         //     successful match for next sibling/expr, it is enough
         //     to continue with matching next exprs.
         if (res)  finalRes = true;

         SX_CHECK (*tmpVisited == *currVisited);

      }

      // --- if no match found then return
      if (finalRes == false)  return false;

      prevMatches = currMatches;
      prevSels    = currSels;

      *currVisited = std::move(*successVisited);

   }
   SxGQuery::SelSet rSet = SxGQuery::SelSet::create ();
   for (auto sIt = prevSels->begin (); sIt != prevSels->end (); ++sIt)  {
      rSet->append ( std::move( *(crossSelections (sels, *sIt)) ) );
   }

   *sels = *rSet;

   *visited = std::move(*currVisited);
   return true;
}

// ---------------------------------------------------------------------------
// single match functions


template<class N,class E,template<class,bool> class GS>
bool SxGQuery::matchInNode (const SxGraph<SxGQPattern>::ConstIterator &it,
                            const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                            const SxGQuery::Selection &sel,
                            const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   ssize_t siz     = sel->getSize ();
   ssize_t visSize = visited->getSize ();
   bool finalRes = false;
   for (ssize_t i = 0; i < gIt.getSizeIn (); ++i)  {
      if(matchRecursive<N,E,GS> (it, gIt.in (i), sel, visited))  {
         // --- atleast one match found
         finalRes = true;
         break;
      }  else  {
         sel->resize (siz);
         visited->resize (visSize);
      }
   }
   return finalRes;
}

template<class N,class E,template<class,bool> class GS>
bool SxGQuery::matchInNodes (const SxGraph<SxGQPattern>::ConstIterator &it,
                             const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                             const SxGQuery::Selection &sel,
                             const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   SxGQuery::Selection resSel  = SxGQuery::Selection::create ();

   SxGQuery::Selection currVisited  = SxGQuery::Selection::create (*visited);
   bool res = true, hasVisited = false;

   for (ssize_t i = 0; i < it.getSizeIn (); ++i)  {
      auto inIt = it.in (i);
      if (!visited->contains (inIt->getId ()))  {
         hasVisited = true;
         res = matchInNode<N,E,GS> (inIt, gIt, resSel, currVisited);
         if (!res)  break;
      }
   }
   if (res && hasVisited)  {
      resSel->append (*sel);
      *sel = *resSel;
      *visited = std::move (*currVisited);
   }
   return res;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchRecursive (const SxGraph<SxGQPattern>::ConstIterator &it,
                          const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                          const SxGQuery::Selection &sel,
                          const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (gIt.isValid ());

   // --- if this node has already been visited
   //     still evaluate/match it so that calling function
   //     can know success/failure. But don't append
   //     it to the result.
   if (visited->contains (it->getId ()))  {
      SxGQuery::Selection tmpSel = SxGQuery::Selection::create ();
      return it->match<N,E,GS> (gIt, tmpSel);
   }
   bool res = false;
   ssize_t siz = sel->getSize ();
   res = it->match<N,E,GS> (gIt, sel);

   if (res)  {
      visited->append (it->getId ());
      res = matchInNodes<N,E,GS> (it, gIt, sel, visited);
   }

   if (res && it.getSizeOut () > 0)  {
      switch (it->getRelType ())  {
         case SxGQPattern::RelationType::OrderedDirect:
            res = matchOutNodesOD<N,E,GS> (it, gIt, sel, visited);
            break;
         case SxGQPattern::RelationType::OrderedIndirect:
            res = matchOutNodesOI<N,E,GS> (it, gIt, sel, visited);
            break;
         case SxGQPattern::RelationType::UnorderedDirect:
            res = matchOutNodesUD<N,E,GS> (it, gIt, sel, visited);
            break;
         case SxGQPattern::RelationType::UnorderedIndirect:
            res = matchOutNodesUI<N,E,GS> (it, gIt, sel, visited);
            break;
      }
   }
   if (!res)  sel->resize (siz);
   return res;
}


template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchOutNodesUDIdx (const SxGraph<SxGQPattern>::ConstIterator &it,
                              const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                              ssize_t childIdx, ssize_t exprChildIdx,
                              const SxGQuery::Selection &sel,
                              const SxGQuery::Selection &visited) const
{
   SX_TRACE ();

   SxGQuery::Selection tmpVisited = SxGQuery::Selection::create (*visited);

   bool res;
   ssize_t siz = sel->getSize ();
   ssize_t visSize = visited->getSize ();
   res = matchRecursive<N,E,GS> (it.out (exprChildIdx), gIt.out (childIdx),
                                 sel, visited);

   if (res == true)  {
      ++exprChildIdx;
      if (exprChildIdx < it.getSizeOut ())  {
         bool res1 = false, res2 = false;
         ssize_t prevSize = 0;
         ssize_t prevVisSize = 0;
         if ((childIdx+1) < gIt.getSizeOut ())  {
            prevSize    = sel->getSize ();
            prevVisSize = visited->getSize ();
            res1 = matchOutNodesUDIdx<N,E,GS> (it, gIt, (childIdx+1),
                                               exprChildIdx, sel, visited);
            if (res1)  {
               return true;
            }  else  {
               SX_CHECK (sel->getSize () == prevSize);
               SX_CHECK (visited->getSize () == prevVisSize);
               SX_UNUSED (prevSize, prevVisSize);
            }
         }

         if ((childIdx-1) >= 0)  {
            prevSize = sel->getSize ();
            prevVisSize = visited->getSize ();
            res2 = matchOutNodesUDIdx<N,E,GS> (it, gIt, (childIdx-1),
                                               exprChildIdx, sel, visited);
            if (res2)  {
               return true;
            }  else  {
               SX_CHECK (sel->getSize () == prevSize);
               SX_CHECK (visited->getSize () == prevVisSize);
               SX_UNUSED (prevSize, prevVisSize);
            }
         }

         res = res1 || res2;
         // --- incase could not match either of next or prev neighbor
         if (!res)  {
            sel->resize (siz);
            visited->resize (visSize);
         }

      }  else  {
         // --- path completed
         return true;
      }
   }

   SX_CHECK (sel->getSize () == siz);
   SX_CHECK (visited->getSize () == visSize);

   return res;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchOutNodesUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                           const SxGQuery::Selection &sel,
                           const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   for (ssize_t i=0; i < gIt.getSizeOut (); ++i)  {
      if (matchOutNodesUDIdx<N,E,GS> (it, gIt, i, 0, sel, visited))
         return true;
   }
   return false;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchOutNodesUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                           const SxGQuery::Selection &sel,
                           const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   ssize_t siz = sel->getSize ();
   ssize_t visSize = sel->getSize ();
   ssize_t outSize = gIt.getSizeOut ();
   ssize_t childIdx = 0;
   ssize_t exprChildIdx = 0;
   bool finalRes = false;

   auto chIt = it.out (exprChildIdx);
   ssize_t prevSize = 0;
   ssize_t prevVisSize = 0;
   for (; childIdx < outSize; ++childIdx)  {
      prevSize    = sel->getSize ();
      prevVisSize = visited->getSize ();
      if (matchRecursive<N,E,GS> (chIt, gIt.out (childIdx), sel, visited))  {
         // --- go to next expr
         exprChildIdx++;
         if (exprChildIdx >= it.getSizeOut ())  {
            finalRes = true;
            break;
         }

         chIt = it.out (exprChildIdx);

         // --- loop will increment it to 0
         childIdx = -1;
      }  else  {
         sel->resize (prevSize);
         visited->resize (prevVisSize);
      }
   }
   if (!finalRes)  {
      sel->resize (siz);
      visited->resize (visSize);
   }

   return finalRes;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchOutNodesODIdx (const SxGraph<SxGQPattern>::ConstIterator &pIt,
                              const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                              const SxGQuery::Selection &sel,
                              ssize_t childIdx,
                              const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   ssize_t siz = sel->getSize ();
   ssize_t visSize = sel->getSize ();
   if ((gIt.getSizeOut ()-childIdx) < pIt.getSizeOut ())  {
      return false;
   }

   for (ssize_t j=0; j < pIt.getSizeOut (); ++j)  {
      auto chIt = pIt.out (j);
      if (!matchRecursive<N,E,GS> (chIt, gIt.out (childIdx), sel, visited))  {
         sel->resize (siz);
         visited->resize (visSize);
         return false;
      }

      ++childIdx;
   }

   return true;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchOutNodesOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                           const SxGQuery::Selection &sel,
                           const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   for (ssize_t i=0; i < gIt.getSizeOut (); ++i)  {
      if (matchOutNodesODIdx<N,E,GS> (it, gIt, sel, i, visited))  return true;
   }

   return false;
}

template<class N,class E,template<class,bool> class GS>
bool
SxGQuery::matchOutNodesOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                           const SxGQuery::Selection &sel,
                           const SxGQuery::Selection &visited) const
{
   SX_TRACE ();
   ssize_t siz = sel->getSize ();
   ssize_t visSize = visited->getSize ();
   ssize_t exprsLeft = it.getSizeOut ();
   ssize_t outSize = gIt.getSizeOut ();

   bool finalRes = false;
   ssize_t i = 0;
   auto chIt = it.out (i);

   ssize_t prevSize = 0;
   ssize_t prevVisSize = 0;
   for (ssize_t j=0; j < gIt.getSizeOut (); ++j)  {

      if ((outSize-j) < exprsLeft)  break;

      prevSize    = sel->getSize ();
      prevVisSize = visited->getSize ();
      if (matchRecursive<N,E,GS> (chIt, gIt.out (j), sel, visited))  {
         ++i;
         if (i >= it.getSizeOut ())  {
            finalRes = true;
            break;
         }

         chIt = it.out (i);
         exprsLeft--;
      }  else  {
         sel->resize (prevSize);
         visited->resize (prevVisSize);
      }
   }

   if (!finalRes)  {
      sel->resize (siz);
      visited->resize (visSize);
   }

   return finalRes;
}
