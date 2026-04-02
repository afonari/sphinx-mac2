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

template<class T>
size_t sxHash (const T &in)
{
   return SxHashFunction::hash (in);
}

// ----------------------------------------------------------------------------
template<class SValue,class SNode,class SContainer,class IT>
SxGraphItState<SValue,SNode,SContainer,IT>::SxGraphItState ()
   : dir(sx::Undefined),
     container(NULL),
     node(NULL),
     selection(),
     idx(-1),
     distMapIdx(-1),
     maxDepth(0),
     pathLen(0),
     distMap(0,64), // where 64 is size of allocation unit
     selectedSizeIn(-1),
     selectedSizeOut(-1)
{
   SX_TRACE ();
}

template<class SValue,class SNode,class SContainer,class IT>
SxGraphItState<SValue,SNode,SContainer,IT>::SxGraphItState (
   sx::Direction dir_, SContainer *c, ssize_t idx_,
   const Selection &selection_)
   : dir(dir_),
     container(c),
     node(NULL),
     selection(selection_),
     idx(idx_),
     distMapIdx(-1),
     maxDepth(0),
     pathLen(0),
     distMap(0,64), // where 64 is size of allocation unit
     selectedSizeIn(-1),
     selectedSizeOut(-1)
{
   SX_TRACE ();
   SX_CHECK (container);

   if (dir != sx::Undefined)  {
      if (idx >= 0)  {
         //distMap << SxPair<ssize_t,ssize_t>(idx, pathLen);
         distMap.append (SxPair<ssize_t,ssize_t>(idx, pathLen));
         visited << idx;
         distMapIdx = -1;
         if (dir == sx::Forward)  procNext ();
         else                     procPrev ();
      }
   }  else  {
      if (selection.getPtr ())  {
         selectionIt = selection->begin ();
         if (selectionIt.isValid ()) ++selectionIt;
      }
      node = (idx >= 0) ? &container->nodes->get (idx) : NULL;
   }
}

template<class SValue,class SNode,class SContainer,class IT>
SxGraphItState<SValue,SNode,SContainer,IT>::SxGraphItState (
     const SxGraphItState &in, sx::ItCopyMode cMode_)
{
   SX_TRACE ();
   copy (in, cMode_);
}

template<class SValue,class SNode,class SContainer,class IT>
SxGraphItState<SValue,SNode,SContainer,IT>::SxGraphItState (
     SxGraphItState &&in, sx::ItCopyMode cMode_) noexcept
{
   SX_TRACE ();
   move (std::move(in), cMode_);
}

template<class SValue,class SNode,class SContainer,class IT>
void SxGraphItState<SValue,SNode,SContainer,IT>::procNext ()
{
   SX_TRACE ();

   if (distMap.getSize () < 1)  {
      idx  = -1;
      node = NULL;
      return;
   }

   distMapIdx += 1;

   if (distMapIdx >= distMap.getSize())  {
      idx = -1;
      node = NULL;
      return;
   }

   idx     = distMap(distMapIdx).key;
   if (idx < 0)  {
      node = NULL;
      return;
   }

   node    = &container->nodes->get (idx);
   pathLen = distMap(distMapIdx).value;
   //distMap.removeFirst ();


   if (dir == sx::Forward)  {
      // --- follow outgoing

      const SxArray<SxPair<ssize_t,ssize_t> > &edges = node->out ();
      ssize_t n = edges.size;
      if (selection.getPtr ())  {
         for (ssize_t i=0; i < n; ++i)  {
            ssize_t outIdx = edges(i).value;
            if (  !visited.contains (outIdx)
               && selection->contains (outIdx))
            {
               //distMap << SxPair<ssize_t,ssize_t>(outIdx,
                 //                                 pathLen + 1);
               distMap.append (SxPair<ssize_t,ssize_t>(outIdx,
                                                       pathLen + 1));
               visited << outIdx;
            }
         }
      }  else  {
         for (ssize_t i=0; i < n; ++i)  {
            ssize_t outIdx = edges(i).value;
            if (!visited.contains (outIdx))  {
               //distMap << SxPair<ssize_t,ssize_t>(outIdx,
                   //                               pathLen + 1);
               distMap.append (SxPair<ssize_t,ssize_t>(outIdx,
                                                       pathLen + 1));
               visited << outIdx;
            }
         }
      }

   }

   if (dir == sx::Backward)  {
      // --- follow incoming

      const SxArray<SxPair<ssize_t,ssize_t> > &edges = node->in ();
      ssize_t n = edges.size;
      if (selection.getPtr ())  {
         for (ssize_t i=0; i < n; ++i)  {
            ssize_t inIdx = edges(i).value;
            if (!visited.contains (inIdx)
                && selection->contains (inIdx))
            {
               //distMap << SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1);
               distMap.append (SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1));
               visited << inIdx;
            }
         }
      }  else  {
         for (ssize_t i=0; i < n; ++i)  {
            ssize_t inIdx = edges(i).value;
            if (!visited.contains (inIdx))  {
               //distMap << SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1);
               distMap.append (SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1));
               visited << inIdx;
            }
         }
      }

   }

}


template<class SValue,class SNode,class SContainer,class IT>
template<class CIT>
void SxGraphItState<SValue,SNode,SContainer,IT>::copy (const CIT &in_,
                                                       sx::ItCopyMode cMode_)
{
   SX_TRACE ();
   container       = in_.container;
   selection       = in_.selection;
   node            = in_.node;
   idx             = in_.idx;
   dir             = in_.dir;

   if (cMode_ == sx::CopyItData)  {
      maxDepth        = 0;
      pathLen         = 0;
      selectedSizeIn  = 0;
      selectedSizeOut = 0;
      distMapIdx      = -1;
      distMap.resize (0);
      visited         = SxSet<ssize_t> ();
      selectedIn      = SxArray<ssize_t> ();
      selectedOut     = SxArray<ssize_t> ();
      selectionIt     = SxList<ssize_t>::ConstIterator ();

      if (dir != sx::Undefined)  {
         if (idx >= 0)  {
            //distMap << SxPair<ssize_t,ssize_t>(idx, pathLen);
            distMap.append (SxPair<ssize_t,ssize_t>(idx, pathLen));
            visited << idx;
            if (dir == sx::Forward)  procNext ();
            else                     procPrev ();
         }
      }  else  {
         if (selection.getPtr ())  {
            selectionIt = selection->begin ();
            if (selectionIt.isValid ()) ++selectionIt;
         }
         node = (idx >= 0) ? &container->nodes->get (idx) : NULL;
      }
   }
   if (cMode_ == sx::CopyItMeta || cMode_ == sx::CopyAll)  {
      maxDepth        = in_.maxDepth;
      pathLen         = in_.pathLen;
      selectedSizeIn  = in_.selectedSizeIn;
      selectedSizeOut = in_.selectedSizeOut;
      distMap         = in_.distMap;
      distMapIdx      = in_.distMapIdx;
      visited         = in_.visited;
      selectedIn      = in_.selectedIn;
      selectedOut     = in_.selectedOut;
      selectionIt     = in_.selectionIt;
   }
}

template<class SValue,class SNode,class SContainer,class IT>
void SxGraphItState<SValue,SNode,SContainer,IT>::move
     (SxGraphItState<SValue,SNode,SContainer,IT> &&in_, sx::ItCopyMode cMode_)
{
   SX_TRACE ();

   container       = in_.container;
   selection       = in_.selection;
   node            = in_.node;
   idx             = in_.idx;
   dir             = in_.dir;

   if (cMode_ == sx::CopyItData)  {
      maxDepth        = 0;
      pathLen         = 0;
      selectedSizeIn  = 0;
      selectedSizeOut = 0;
      distMapIdx      = -1;
      distMap.resize (0);
      visited         = SxSet<ssize_t> ();
      selectedIn      = SxArray<ssize_t> ();
      selectedOut     = SxArray<ssize_t> ();
      selectionIt     = SxList<ssize_t>::ConstIterator ();

      if (dir != sx::Undefined)  {
         if (idx >= 0)  {
            //distMap << SxPair<ssize_t,ssize_t>(idx, pathLen);
            distMap.append (SxPair<ssize_t,ssize_t>(idx, pathLen));
            visited << idx;
            if (dir == sx::Forward)  procNext ();
            else                     procPrev ();
         }
      }  else  {
         if (selection.getPtr ())  {
            selectionIt = selection->begin ();
            if (selectionIt.isValid ()) ++selectionIt;
         }
         node = (idx >= 0) ? &container->nodes->get (idx) : NULL;
      }
   }

   if (cMode_ == sx::CopyItMeta || cMode_ == sx::CopyAll)  {
      maxDepth        = in_.maxDepth;
      pathLen         = in_.pathLen;
      selectedSizeIn  = in_.selectedSizeIn;
      selectedSizeOut = in_.selectedSizeOut;
      distMap         = std::move (in_.distMap);
      distMapIdx      = in_.distMapIdx;
      visited         = std::move (in_.visited);
      selectedIn      = std::move (in_.selectedIn);
      selectedOut     = std::move (in_.selectedOut);
      selectionIt     = std::move (in_.selectionIt);
   }

   in_.dir             = sx::Undefined;
   in_.container       = NULL;
   in_.node            = NULL;
   in_.selection       = Selection ();// copy selection
   in_.idx             = -1;
   in_.maxDepth        = 0;
   in_.pathLen         = 0;
   in_.selectedSizeIn  = -1;
   in_.selectedSizeOut = -1;

}

template<class SValue,class SNode,class SContainer,class IT>
IT SxGraphItState<SValue,SNode,SContainer,IT>::insertElem (ssize_t newPos,
                                                           const SValue &elem)
{
   SX_TRACE ();
   SX_CHECK (container);
   SX_UNUSED (newPos);
   return container->createNode (elem);
}

template<class SValue,class SNode,class SContainer,class IT>
IT SxGraphItState<SValue,SNode,SContainer,IT>::appendElem (const SValue &elem)
{
   SX_TRACE ();
   SX_CHECK (container);
   return container->createNode (elem);
}

template<class SValue,class SNode,class SContainer,class IT>
IT SxGraphItState<SValue,SNode,SContainer,IT>::prependElem (const SValue &elem)
{
   SX_TRACE ();
   SX_CHECK (container);
   return container->createNode (elem);
}

template<class SValue,class SNode,class SContainer,class IT>
void SxGraphItState<SValue,SNode,SContainer,IT>::next ()
{
   SX_CHECK (container);

   if (dir == sx::Backward)  {
      if (distMap.getSize () < 1)  {
         idx  = -1;
         node = NULL;
         return;
      }

      distMapIdx -= 1;

      if (distMapIdx < 0)  {
         idx = -1;
         node = NULL;
         return;
      }

      idx     = distMap(distMapIdx).key;
      if (idx < 0)  {
         node = NULL;
         return;
      }
      node    = &container->nodes->get (idx);
      pathLen = distMap(distMapIdx).value;
   }  else if (dir != sx::Undefined)  {
      procNext ();
   }  else if (selection.getPtr ())  {
      if (selectionIt.isValid ())  {
         idx = *selectionIt;
         node = &container->nodes->get (idx);
         ++selectionIt;
      }  else  {
         idx = -1;
         node = NULL;
      }
   }  else  {
      //idx = node->next;
      idx = container->nodes->next (idx);
      node = (idx >= 0) ? &container->nodes->get (idx) : NULL;
   }
   // --- reset neighbors in selection
   selectedSizeIn = -1;
   selectedSizeOut = -1;

   if (maxDepth > 0 && getDepth() > maxDepth)  {
      node = NULL;
      idx = -1;
   }
}

template<class SValue,class SNode,class SContainer,class IT>
void SxGraphItState<SValue,SNode,SContainer,IT>::procPrev ()
{
   SX_TRACE ();

   if (distMap.getSize () < 1)  {
      idx  = -1;
      node = NULL;
      return;
   }

   distMapIdx += 1;

   if (distMapIdx >= distMap.getSize())  {
      idx = -1;
      node = NULL;
      return;
   }

   idx     = distMap(distMapIdx).key;
   if (idx < 0)  {
      node = NULL;
      return;
   }

   node    = &container->nodes->get (idx);
   pathLen = distMap(distMapIdx).value;

   if (dir == sx::Backward)  {
      // --- follow incoming

      const SxArray<SxPair<ssize_t,ssize_t> > &edges = node->in ();
      ssize_t n = edges.size;
      if (selection.getPtr ())  {
         for (ssize_t i=0; i < n; ++i)  {
            ssize_t inIdx = edges(i).value;
            if (!visited.contains (inIdx)
                && selection->contains (inIdx))
            {
               //distMap << SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1);
               distMap.append (SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1));
               visited << inIdx;
            }
         }
      }  else  {
         for (ssize_t i=0; i < n; ++i)  {
            ssize_t inIdx = edges(i).value;
            if (!visited.contains (inIdx))  {
               //distMap << SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1);
               distMap.append (SxPair<ssize_t,ssize_t>(inIdx, pathLen + 1));
               visited << inIdx;
            }
         }
      }

   }

}


template<class SValue,class SNode,class SContainer,class IT>
void SxGraphItState<SValue,SNode,SContainer,IT>::prev ()
{
   SX_CHECK (container);
   //SX_CHECK (dir != sx::Forward);
   if (dir == sx::Forward)  {
      if (distMap.getSize () < 1)  {
         idx  = -1;
         node = NULL;
         return;
      }

      distMapIdx -= 1;

      if (distMapIdx < 0)  {
         idx = -1;
         node = NULL;
         return;
      }

      idx     = distMap(distMapIdx).key;
      if (idx < 0)  {
         node = NULL;
         return;
      }
      node    = &container->nodes->get (idx);
      pathLen = distMap(distMapIdx).value;
   }  else if (dir != sx::Undefined)  {
      procPrev ();
   }  else if (selection.getPtr ())  {
      if (selectionIt.isValid ())  {
         idx = *selectionIt;
         node = &container->nodes->get (idx);
         ++selectionIt;
      }  else  {
         idx = -1;
         node = NULL;
      }
   }  else  {
      idx = container->nodes->prev (idx);
      node = (idx >= 0) ? &container->nodes->get (idx) : NULL;
   }
   // --- reset neighbors in selection
   selectedSizeIn = -1;
   selectedSizeOut = -1;

   if (maxDepth > 0 && getDepth() > maxDepth)  {
      node = NULL;
      idx = -1;
   }

}

template<class SValue,class SNode,class SContainer,class IT>
SValue &SxGraphItState<SValue,SNode,SContainer,IT>::getRef ()
{
   SX_TRACE ();
   SX_CHECK (node);
   return node->element;
}

template<class SValue,class SNode,class SContainer,class IT>
SValue *SxGraphItState<SValue,SNode,SContainer,IT>::getPtr ()
{
   SX_TRACE ();
   SX_CHECK (node);
   return &node->element;
}

template<class SValue,class SNode,class SContainer,class IT>
bool SxGraphItState<SValue,SNode,SContainer,IT>::equal (const IT &it) const
{
   SX_TRACE ();
   SX_CHECK (container == it.container || !container || !it.container);
   return idx == it.idx;
}

template<class SValue,class SNode,class SContainer,class IT>
typename SxGraphItState<SValue,SNode,SContainer,IT>::Selection
SxGraphItState<SValue,SNode,SContainer,IT>::getUnvisited () const
{
   SX_TRACE ();
   SX_CHECK (container);
   Selection result;

   IT it = container->begin ();
   while (it.isValid ())  {
      ssize_t nodeId = it.idx;
      if (!visited.contains (nodeId))  {
         if (!result.getPtr ())
            result = Selection::create ();
         result->append (nodeId);
      }
      ++it;
   }
   return result;
}


template<class SValue,class SNode,class SContainer,class IT>
IT SxGraphItState<SValue,SNode,SContainer,IT>::in (ssize_t idx_) const
{
   SX_TRACE ();
   SX_CHECK (node);

   if (selection.getPtr ())  {
      if (selectedSizeIn < 0)  {
         // --- build cache of selected neighbors
         selectedSizeIn = 0;
         const SxArray<SxPair<ssize_t,ssize_t> > &edgesIn = node->in ();
         ssize_t n = edgesIn.size;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesIn(i).value))  {
               selectedSizeIn++;
            }
         }
         selectedIn.resize (selectedSizeIn);
         selectedSizeIn = 0;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesIn(i).value))  {
               selectedIn(selectedSizeIn) = i;
               selectedSizeIn++;
            }
         }
      }
      return IT (dir, container,
                 node->in ()(selectedIn(idx_)).value,
                 selection);
   }
   return IT (dir, container, node->in ()(idx_).value,
              selection);

}

template<class SValue,class SNode,class SContainer,class IT>
IT SxGraphItState<SValue,SNode,SContainer,IT>::out (ssize_t idx_) const
{
   SX_TRACE ();
   SX_CHECK (node);

   if (selection.getPtr ())  {
      if (selectedSizeOut < 0)  {
         // --- build cache of selected neighbors
         selectedSizeOut = 0;
         const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = node->out ();
         ssize_t n = edgesOut.size;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesOut(i).value))  {
               selectedSizeOut++;
            }
         }
         selectedOut.resize (selectedSizeOut);
         selectedSizeOut = 0;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesOut(i).value))  {
               selectedOut(selectedSizeOut) = i;
               selectedSizeOut++;
            }
         }
      }
      return IT (dir, container,
                 node->out ()(selectedOut(idx_)).value,
                 selection);
   }
   return IT(dir, container, node->out ()(idx_).value,
             selection);

}

template<class SValue,class SNode,class SContainer,class IT>
IT SxGraphItState<SValue,SNode,SContainer,IT>::neighbors (ssize_t maxDepth_) const
{
   SX_TRACE ();
   SX_CHECK (maxDepth_ >= 0, maxDepth_);
   IT res; res.copy (*this, sx::CopyItData); // parent
   res.maxDepth = maxDepth_;
   // iterate to first child
   if (res.getDirection () == sx::Forward)  ++res;
   else                                     --res;
   return res;
}

template<class SValue,class SNode,class SContainer,class IT>
IT SxGraphItState<SValue,SNode,SContainer,IT>::hops (ssize_t maxDepth_) const
{
   SX_TRACE ();
   SX_CHECK (maxDepth_ >= 0, maxDepth_);
   IT res; res.copy (*this, sx::CopyItData); // parent
   res.maxDepth = maxDepth_;
   return res;
}

template<class SValue,class SNode,class SContainer,class IT>
SValue &SxGraphItState<SValue,SNode,SContainer,IT>::getIn (ssize_t idx_)
{
   SX_TRACE ();
   SX_CHECK (node);

   if (selection.getPtr ())  {
      if (selectedSizeIn < 0)  {
         // --- build cache of selected neighbors
         selectedSizeIn = 0;
         const SxArray<SxPair<ssize_t,ssize_t> > &edgesIn = node->in ();
         ssize_t n = edgesIn.size;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesIn(i).value))  {
               selectedSizeIn++;
            }
         }
         selectedIn.resize (selectedSizeIn);
         selectedSizeIn = 0;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesIn(i).value))  {
               selectedIn(selectedSizeIn) = i;
               selectedSizeIn++;
            }
         }
      }
      return container->nodes->get (node->out ()(
                                    selectedIn(idx_)).value).element;
   }
   return container->nodes->get (node->in ()(idx_).value).element;

}

template<class SValue,class SNode,class SContainer,class IT>
SValue &SxGraphItState<SValue,SNode,SContainer,IT>::getOut (ssize_t idx_)
{
   SX_TRACE ();
   SX_CHECK (node);

   if (selection.getPtr ())  {
      if (selectedSizeOut < 0)  {
         // --- build cache of selected neighbors
         selectedSizeOut = 0;
         const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = node->out ();
         ssize_t n = edgesOut.size;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesOut(i).value))  {
               selectedSizeOut++;
            }
         }
         selectedOut.resize (selectedSizeOut);
         selectedSizeOut = 0;
         for (ssize_t i=0; i < n; ++i)  {
            if (selection->contains (edgesOut(i).value))  {
               selectedOut(selectedSizeOut) = i;
               selectedSizeOut++;
            }
         }
      }
      return container->nodes->get (node->out ()(
                                    selectedOut(idx_)).value).element;
   }
   return container->nodes->get (node->out ()(idx_).value).element;

}

template<class SValue,class SNode,class SContainer,class IT>
ssize_t SxGraphItState<SValue,SNode,SContainer,IT>::getSizeIn () const
{
   SX_TRACE ();
   SX_CHECK (node);
   if (selection.getPtr ())  {
      // --- count selected neighbors
      ssize_t nSelected = 0;

      const SxArray<SxPair<ssize_t,ssize_t> > &edgesIn = node->in ();
      ssize_t n = edgesIn.size;
      for (ssize_t i=0; i < n; ++i)  {
         if (selection->contains (edgesIn(i).value))  {
            nSelected++;
         }
      }

      return nSelected;
   }
   return node->in ().size;
}

template<class SValue,class SNode,class SContainer,class IT>
ssize_t SxGraphItState<SValue,SNode,SContainer,IT>::getSizeOut () const
{
   SX_TRACE ();
   SX_CHECK (node);
   if (selection.getPtr ())  {
      // --- count selected neighbors
      ssize_t nSelected = 0;

      const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = node->out ();
      ssize_t n = edgesOut.size;
      for (ssize_t i=0; i < n; ++i)  {
         if (selection->contains (edgesOut(i).value))  {
            nSelected++;
         }
      }

      return nSelected;
   }
   return node->out ().size;
}

template<class SValue,class SNode,class SContainer,class IT>
ssize_t SxGraphItState<SValue,SNode,SContainer,IT>::findEdgeOut (ssize_t idxTo)
{
   SX_TRACE ();
   SX_CHECK (idx >= 0);
   return container->findEdgeOut (idx, idxTo);
}

template<class SValue,class SNode,class SContainer,class IT>
ssize_t SxGraphItState<SValue,SNode,SContainer,IT>::findEdgeIn (ssize_t idxTo)
{
   SX_TRACE ();
   SX_CHECK (idx >= 0);
   return container->findEdgeIn (idx, idxTo);
}

// ----------------------------------------------------------------------------

template<class T,bool isNS>
SxGraphStorage<T,isNS>::SxGraphStorage ()
   : elems(0, 64),
     firstElement(-1),
     lastElement(-1),
     firstFree(-1),
     nElems(0),
     hashSize(0),
     hashUsed(0),
     hashBound(0)
{
   SX_TRACE ();
}

template<class T,bool isNS>
SxGraphStorage<T,isNS>::SxGraphStorage (const SxGraphStorage<T,isNS> &in)
   : elems(in.elems),
     hashTable(in.hashTable),
     firstElement(in.firstElement),
     lastElement(in.lastElement),
     firstFree(in.firstFree),
     nElems(in.nElems),
     hashSize(in.hashSize),
     hashUsed(in.hashUsed),
     hashBound(in.hashBound)
{
   SX_TRACE ();
}

template<class T,bool isNS>
SxGraphStorage<T,isNS>::SxGraphStorage (SxGraphStorage<T,isNS> &&in)
   : elems(std::move(in.elems)),
     hashTable(std::move(in.hashTable)),
     firstElement(in.firstElement),
     lastElement(in.lastElement),
     firstFree(in.firstFree),
     nElems(in.nElems),
     hashSize(in.hashSize),
     hashUsed(in.hashUsed),
     hashBound(in.hashBound)
{
   SX_TRACE ();
   in.firstElement = -1;
   in.lastElement = -1;
   in.firstFree = -1;
   in.nElems = 0;
   in.hashSize = 0;
   in.hashUsed = 0;
   in.hashBound = 0;
}

template<class T,bool isNS>
SxGraphStorage<T,isNS> &
SxGraphStorage<T,isNS>::operator= (const SxGraphStorage<T,isNS> &in)
{
   SX_TRACE ();
   if (this == &in)  return *this;

   elems        = in.elems;
   hashTable    = in.hashTable;
   firstElement = in.firstElement;
   lastElement  = in.lastElement;
   firstFree    = in.firstFree;
   nElems       = in.nElems;

   hashSize     = in.hashSize;
   hashUsed     = in.hashUsed;
   hashBound    = in.hashBound;

   return *this;
}

template<class T,bool isNS>
SxGraphStorage<T,isNS> &
SxGraphStorage<T,isNS>::operator= (SxGraphStorage<T,isNS> &&in)
{
   SX_TRACE ();
   if (this == &in)  return *this;

   elems        = std::move(in.elems);
   hashTable    = std::move(in.hashTable);
   firstElement = in.firstElement;
   lastElement  = in.lastElement;
   firstFree    = in.firstFree;
   nElems       = in.nElems;

   hashSize     = in.hashSize;
   hashUsed     = in.hashUsed;
   hashBound    = in.hashBound;

   in.firstElement = -1;
   in.lastElement  = -1;
   in.firstFree    = -1;
   in.nElems       = 0;
   in.hashSize     = 0;
   in.hashUsed     = 0;
   in.hashBound    = 0;

   return *this;
}

template<class T,bool isNS>
SxGraphStorage<T,isNS>::~SxGraphStorage ()
{
   removeAll ();
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::addElem (const T &elem)
{
   SX_TRACE ();
   ssize_t idxHash = 0;
   if (!hashContains (elem, &idxHash))  {
      ssize_t idx = 0;
      if (firstFree >= 0)  {
         idx = firstFree;
         firstFree = elems(idx).next;
      }  else  {
         //idx = static_cast<ssize_t>(elems.getSize ());
         idx = static_cast<ssize_t>(elems.size);
         elems.resize (elems.size + 1);
      }
      Elem &entry   = elems(idx);
      entry.element = elem;
      entry.prev    = lastElement;
      entry.next    = -1;

      if (lastElement >= 0)  elems(lastElement).next = idx;
      lastElement = idx;
      if (firstElement < 0)  firstElement = idx;
      nElems++;

      hashTable.elements[idxHash] = lastElement;
      return idx;
   }
   return -1;
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::addElem (T &&elem)
{
   SX_TRACE ();
   ssize_t idxHash = 0;

   if (!hashContains (elem, &idxHash))  {
      ssize_t idx = 0;
      if (firstFree >= 0)  {
         idx = firstFree;
         firstFree = elems(idx).next;
      }  else  {
         //idx = static_cast<ssize_t>(elems.getSize ());
         idx = static_cast<ssize_t>(elems.size);
         elems.resize (elems.size + 1);
      }
      Elem &entry   = elems(idx);
      entry.element = std::move(elem);
      entry.prev    = lastElement;
      entry.next    = -1;

      if (lastElement >= 0)  elems(lastElement).next = idx;
      lastElement = idx;
      if (firstElement < 0)  firstElement = idx;
      nElems++;

      hashTable.elements[idxHash] = lastElement;
      return idx;
   }
   return -1;
}

template<class T,bool isNS>
void SxGraphStorage<T,isNS>::removeElem (ssize_t idx)
{
   SX_TRACE ();
   //SX_CHECK (elems.containsKey (idx));
   if (containsElem (idx))  {
      ssize_t prevPtr = elems(idx).prev;
      ssize_t nextPtr = elems(idx).next;
      ssize_t hashIdx = hashFindPos (elems(idx).element);
      if (hashIdx >= 0)  hashTable.elements[hashIdx] = -2; // deleted
      elems(idx) = Elem(); // destroy

      if (prevPtr >= 0) elems(prevPtr).next = nextPtr;
      else              firstElement = nextPtr;
      if (nextPtr >= 0) elems(nextPtr).prev = prevPtr;
      else              lastElement = prevPtr;

      elems(idx).next = firstFree;
      elems(idx).prev = -1;
      firstFree = idx;
      nElems -= 1;

   }
}

template<class T,bool isNS>
bool SxGraphStorage<T,isNS>::containsElem (ssize_t idx) const
{
   SX_TRACE ();
   if (idx >=0 && idx < elems.size)  {
      ssize_t freeIdx = firstFree;
      while (freeIdx >= 0)  {
         if (freeIdx == idx)  {
            return false;
         }
         freeIdx = elems(freeIdx).next;
      }
      return true;
   }
   return false;
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::getIdx (const T &elem)
{
   SX_TRACE ();
   ssize_t idx = hashFindPos (elem);
   return (idx >= 0) ? hashTable.elements[idx] : -1;
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::getFirstIdx () const
{
   SX_TRACE ();
   return firstElement;
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::getNextIdx (ssize_t currIdx) const
{
   SX_TRACE ();
   SX_CHECK (currIdx >= 0);
   SX_CHECK (containsElem (currIdx));
   return elems(currIdx).next;
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::getPrevIdx (ssize_t currIdx) const
{
   SX_TRACE ();
   SX_CHECK (currIdx >= 0);
   SX_CHECK (containsElem (currIdx));
   return elems(currIdx).prev;
}

template<class T,bool isNS>
typename SxGraphStorage<T,isNS>::Elem &
SxGraphStorage<T,isNS>::getElem (ssize_t idx)
{
   SX_TRACE ();
   SX_CHECK (containsElem (idx));
   return  elems(idx);
}

template<class T,bool isNS>
const typename SxGraphStorage<T,isNS>::Elem &
SxGraphStorage<T,isNS>::getElem (ssize_t idx) const
{
   SX_TRACE ();
   SX_CHECK (containsElem (idx));
   return  elems(idx);
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::hashFindPos (const T &elem) const
{
   SX_TRACE ();
   if (hashSize > 0)  {
      size_t h = sxHash (elem);
      size_t tblSize = static_cast<size_t>(hashTable.size);
      size_t idx = h % tblSize;
      size_t step = h % (tblSize - 1) + 1;
      size_t lastIdx = idx;

      do {
         if (hashTable.elements[idx] == -1)  {
            return -1;
         }  else if (hashTable.elements[idx] != -2 &&
                     elems((ssize_t)hashTable.elements[idx]).element == elem)
         {
            return static_cast<ssize_t>(idx);
         }
         idx = (idx + step) % tblSize;
      } while (idx != lastIdx);
   }
   return -1;
}

template<class T,bool isNS>
bool SxGraphStorage<T,isNS>::hashContains (const T &elem_, ssize_t *result_)
{
   SX_TRACE ();
   SX_CHECK (result_);

   if (hashUsed >= hashBound)  {
      if (hashSize * 2 < hashTable.size)  {
         hashResize (hashTable.size - 1);
      }  else  {
         hashResize (hashTable.size + 1);
      }
   }

   size_t h = sxHash (elem_);
   size_t tblSize = static_cast<size_t>(hashTable.size);
   size_t idx = h % tblSize;
   size_t step = (h % (tblSize - 1)) + 1;
   size_t lastIdx = idx;
   size_t replaceIdx = tblSize;

   do {
      if (hashTable.elements[idx] == -1)  {
         // --- insert new element
         if (replaceIdx != tblSize)  {
            idx = replaceIdx;
         }  else  {
            hashUsed++;
         }
         hashSize++;
         *result_ = static_cast<ssize_t>(idx);
         return false;
      }  else if (hashTable.elements[idx] == -2)  {
         replaceIdx = idx;
      }  else if (elems((ssize_t)hashTable.elements[idx]).element == elem_)  {
         *result_ = static_cast<ssize_t>(idx);
         return true;
      }
      idx = (idx + step) % tblSize;
   } while (idx != lastIdx);

   throw SxException ("Can not insert.", __FILE__, __LINE__);
}

template<class T,bool isNS>
void SxGraphStorage<T,isNS>::hashResize (ssize_t newSize_)
{
   SX_TRACE ();
   if (newSize_ >= 1610612741)  {
      throw SxException ("HashTable hashSize is too big", __FILE__, __LINE__);
   }

   // --- http://planetmath.org/encyclopedia/GoodHashTablePrimes.html
   ssize_t prime[30] = {0, 3, 11, 23, 53, 97, 193, 389, 769, 1543, 3079, 6151,
                        12289, 24593, 49157, 98317, 196613, 393241, 786433,
                        1572869, 3145739, 6291469, 12582917, 25165843, 50331653,
                        100663319, 201326611, 402653189, 805306457, 1610612741};
   size_t idx = 0;
   while (prime[idx] <= newSize_)  {
      idx++;
   }
   newSize_ = prime[idx];
   ssize_t newBound = (ssize_t)((double)newSize_ * 0.77 + 0.5);
   ssize_t n = hashTable.size;
   ssize_t i;

   if (hashSize < newBound)  {
      SxArray<ssize_t> newTable(newSize_);
      if (newSize_ > 0)  {
         newTable.set (-1);
      }

      for (i = 0; i < n; i++)  {
         if (hashTable.elements[i] >= 0)  {
            size_t h = sxHash (elems(hashTable.elements[i]).element);
            size_t newSiz_ = static_cast<size_t>(newSize_);
            idx = (h % newSiz_);
            size_t step = (h % (newSiz_ - 1)) + 1;
            size_t lastIdx = idx;

            while (newTable.elements[idx] >= 0)  {
               idx = (idx + step) % newSiz_;
               if (idx == lastIdx)  {
                  SX_EXIT;
               }
            }
            newTable.elements[idx] = hashTable.elements[i];
         }
      }

      hashUsed = hashSize;
      hashBound = newBound;
      hashTable.resize (newSize_, false);
      size_t len = sizeof(ssize_t) * (size_t)hashTable.size;
      ::memcpy (hashTable.elements, newTable.elements, len);
   }
}

template<class T,bool isNS>
ssize_t SxGraphStorage<T,isNS>::getNElems () const
{
   SX_TRACE ();
   return nElems;
}

template<class T,bool isNS>
size_t SxGraphStorage<T,isNS>::getSizeBytes () const
{
   SX_TRACE ();
   SX_EXIT;
   return 0;
}

template<class T,bool isNS>
void SxGraphStorage<T,isNS>::removeElems ()
{
   SX_TRACE ();
   elems.removeAll ();

   firstElement = -1;
   lastElement = -1;
   firstFree = -1;
   nElems = 0;

   hashTable.removeAll ();
   hashSize = 0;
   hashUsed = 0;
   hashBound = 0;
}

// ----------------------------------------------------------------------------
template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxGraph<N,E,GS,ItPair>::SxGraph ()
   : SxThis<SxGraph<N,E,GS,ItPair> > (),
     nodes(SxPtr<GS<N,true> >::create()),
     edges(SxPtr<GS<E,false> >::create())
{
   SX_TRACE ();
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxGraph<N,E,GS,ItPair>::SxGraph (const SxGraph<N,E,GS,ItPair> &in)
   : SxThis<SxGraph<N,E,GS,ItPair> > (),
     nodes(in.nodes),
     edges(in.edges)
{
   SX_TRACE ();
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxGraph<N,E,GS,ItPair>::SxGraph (SxGraph<N,E,GS,ItPair> &&in)
   : SxThis<SxGraph<N,E,GS,ItPair> > (),
     nodes(in.nodes),
     edges(in.edges)
{
   SX_TRACE ();
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxGraph<N,E,GS,ItPair>::~SxGraph ()
{
   SX_TRACE ();
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxGraph<N,E,GS,ItPair> &
SxGraph<N,E,GS,ItPair>::operator= (const SxGraph<N,E,GS,ItPair> &in)
{
   SX_TRACE ();
   if (this == &in)  return *this;

   nodes = in.nodes;
   edges = in.edges;

   return *this;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxGraph<N,E,GS,ItPair> &
SxGraph<N,E,GS,ItPair>::operator= (SxGraph<N,E,GS,ItPair> &&in) noexcept
{
   SX_TRACE ();
   if (this == &in)  return *this;

   nodes = in.nodes;
   edges = in.edges;

   return *this;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::createNode (const N &elem)
{
   SX_TRACE ();

   nodes->add (elem);
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator (sx::Forward, this,
                                                          findPos (elem),
                                                          Selection());
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::createNode (N &&elem)
{
   SX_TRACE ();

   nodes->add (std::move(elem));
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator (sx::Forward, this,
                                                          findPos (elem),
                                                          Selection());
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::createEdge (const N &from, const N &to,
                                         const E &elem)
{
   SX_TRACE ();
   using namespace sx;

   ssize_t idxFrom = findPos (from);
   SX_CHECK (idxFrom >= 0, idxFrom);

   ssize_t idxTo   = findPos (to);
   SX_CHECK (idxTo >= 0, idxTo);

   if (hasOutEdge (idxFrom, idxTo)) return;

   //handleEdgeEvent (BeforeEdgeInEvent | BeforeEdgeOutEvent, idxFrom, idxTo);

   ssize_t idx = edges->add (elem);
   if (idx >= 0)  {
      ET    &e = edges->get (idx);
      e.in  () = idxFrom;
      e.out () = idxTo;

      nodes->get (idxFrom).out ().append (SxPair<ssize_t,ssize_t>(idx, idxTo));
      nodes->get (idxTo).in ().append (SxPair<ssize_t,ssize_t>(idx, idxFrom));

      //handleEdgeEvent (sx::AfterEdgeInEvent | sx::AfterEdgeOutEvent,
        //               idxFrom, idxTo);
   }
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::createEdge (typename SxGraph<N,E,GS,ItPair
                                                       >::ConstIterator &fromIt,
                                         typename SxGraph<N,E,GS,ItPair
                                                         >::ConstIterator &toIt,
                                         const E &elem)
{
   SX_TRACE ();
   using namespace sx;

   ssize_t idxFrom = fromIt.getIdx ();
   SX_CHECK (idxFrom >= 0, idxFrom);

   ssize_t idxTo   = toIt.getIdx ();
   SX_CHECK (idxTo >= 0, idxTo);

   if (hasOutEdge (idxFrom, idxTo)) return;

   //handleEdgeEvent (BeforeEdgeInEvent | BeforeEdgeOutEvent, idxFrom, idxTo);

   ssize_t idx = edges->add (elem);
   if (idx >= 0)  {
      ET    &e = edges->get (idx);
      e.in  () = idxFrom;
      e.out () = idxTo;

      nodes->get (idxFrom).out ().append (SxPair<ssize_t,ssize_t>(idx, idxTo));
      nodes->get (idxTo).in ().append (SxPair<ssize_t,ssize_t>(idx, idxFrom));

      //handleEdgeEvent (sx::AfterEdgeInEvent | sx::AfterEdgeOutEvent,
        //               idxFrom, idxTo);
   }
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::removeNode (ssize_t nodeIdx)
{
   SX_TRACE ();
   // --- remove edges

   NT &node = nodes->get (nodeIdx);
   const SxArray<SxPair<ssize_t,ssize_t> > edgesOut = node.out ();
   for (ssize_t i = 0; i < edgesOut.size; ++i)  {
      unlinkIn (nodeIdx, edgesOut.elements[i].value); // in-->out
   }
   const SxArray<SxPair<ssize_t, ssize_t> > edgesIn = node.in ();
   for (ssize_t i = 0; i < edgesIn.size; ++i)  {
      unlinkOut (edgesIn.elements[i].value, nodeIdx);
   }

   // --- remove node
   nodes->remove (nodeIdx);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::removeElement (const N &elem)
{
   SX_TRACE ();
   ssize_t ptr = findPos (elem);
   SX_CHECK (ptr >= 0, ptr);
   removeNode (ptr);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::removeEdge (const N &from, const N &to)
{
   SX_TRACE ();
   ssize_t idxFrom = findPos (from);
   SX_CHECK (idxFrom >= 0, idxFrom);
   ssize_t idxTo   = findPos (to);
   SX_CHECK (idxTo >= 0, idxTo);

   SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get (idxFrom).out ();
   SxList<ssize_t> edgeIdx;
   for (ssize_t i = 0; i < edgesOut.size; ++i)  {
      if (idxTo == edgesOut.elements[i].value)  {
         edgeIdx.append (edgesOut.elements[i].key);
      }
   }

   unlinkOut (idxFrom, idxTo);
   unlinkIn  (idxFrom, idxTo);

   for (auto it = edgeIdx.begin (); it != edgeIdx.end (); ++it)  {
      edges->remove (*it);
   }
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
bool SxGraph<N,E,GS,ItPair>::hasOutEdge (ssize_t idxFrom,
                                         ssize_t idxTo) const
{
   SX_TRACE ();
   return findEdgeOut (idxFrom, idxTo) >= 0;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
ssize_t SxGraph<N,E,GS,ItPair>::findEdgeOut (ssize_t idxFrom,
                                             ssize_t idxTo) const
{
   SX_TRACE ();

   const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get (idxFrom).out ();
   const SxArray<SxPair<ssize_t,ssize_t> > &edgesIn  = nodes->get (idxTo).in ();

   if (edgesOut.size < edgesIn.size)  {
      for (ssize_t i = 0; i < edgesOut.size; ++i)  {
         if (edgesOut.elements[i].value == idxTo)  {
            return edgesOut.elements[i].key;
         }
      }
   }  else  {
      for (ssize_t i = 0; i < (edgesIn.size); ++i)  {
         if (edgesIn.elements[i].value == idxFrom)
            return edgesIn.elements[i].key;
      }
   }

   return -1;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
ssize_t SxGraph<N,E,GS,ItPair>::findEdgeIn (ssize_t idxFrom,
                                            ssize_t idxTo) const
{
   SX_TRACE ();

   const SxArray<SxPair<ssize_t,ssize_t> > &edgesIn  = nodes->get(idxFrom).in ();
   const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get(idxTo).out ();

   if (edgesIn.size < edgesOut.size)  {
      for (ssize_t i = 0; i < (edgesIn.size); ++i)  {
         if (edgesIn.elements[i].value == idxTo)
            return edgesIn.elements[i].key;
      }
   }  else  {
      for (ssize_t i = 0; i < edgesOut.size; ++i)  {
         if (edgesOut.elements[i].value == idxFrom)  {
            return edgesOut.elements[i].key;
         }
      }
   }

   return -1;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::unlinkOut (ssize_t idxFrom, ssize_t idxTo)
{
   SX_TRACE ();
   // --- remove all outgoing
   SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get(idxFrom).out ();
   ssize_t removeNEdges = 0;
   ssize_t idx = -1;
   const ssize_t n = edgesOut.size;
   for (ssize_t i=0; i < n; ++i)  {
      if (idxTo != edgesOut.elements[i].value)  {
         if (idx >= 0)  {
            // --- shift to empty place from removed edges
            edgesOut.elements[idx] = edgesOut.elements[i];
            idx++;
         }
      }  else  {
         removeNEdges++;
         if (idx < 0) idx = i; // insert place
      }
   }
   if (removeNEdges > 0)  {
      edgesOut.remove (n - removeNEdges, removeNEdges);
   }

}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::unlinkIn (ssize_t idxFrom, ssize_t idxTo)
{
   SX_TRACE ();
   // --- remove incomming
   SxArray<SxPair<ssize_t,ssize_t> > &edgesIn = nodes->get (idxTo).in ();

   ssize_t removeNEdges = 0;
   ssize_t idx = -1;
   const ssize_t n = edgesIn.size;
   for (ssize_t i=0; i < n; ++i)  {
      if (idxFrom != edgesIn.elements[i].value)  {
         if (idx >= 0)  {
            // --- shift to empty place from removed edges
            edgesIn.elements[idx] = edgesIn.elements[i];
            idx++;
         }
      }  else  {
         removeNEdges++;
         if (idx < 0) idx = i; // insert place
      }
   }
   if (removeNEdges > 0)  {
      edgesIn.remove (n - removeNEdges, removeNEdges);
   }

}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::removeSelection (const Selection &selection)
{
   SX_TRACE ();
   if (selection.getPtr ())  {
      SxList<ssize_t>::ConstIterator it = selection->begin ();
      while (it.isValid ())  {
         // --- invalid selection
         //     (*it) < 0 || (*it) >= nodes.getSize ()
         //     || (*it) is in firstFree list
         removeNode (*it);
         ++it;
      }
   }
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::removeAll ()
{
   SX_TRACE ();
   nodes->removeAll ();
   edges->removeAll ();
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
bool SxGraph<N,E,GS,ItPair>::containsNode (const N &elem) const
{
   SX_TRACE ();
   return findPos (elem) >= 0;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
ssize_t SxGraph<N,E,GS,ItPair>::findPath (const N &from, const N &to,
                                          ssize_t maxPathLen,
                                          const Selection &selection) const
{
   SX_TRACE ();
   ssize_t idxFrom = findPos (from);
   if (idxFrom < 0 || (selection.getPtr () && !selection->contains (idxFrom)))
      return 0;
   ssize_t idxTo = findPos (to);
   if (idxTo < 0 || (selection.getPtr () && !selection->contains (idxTo)))
      return 0;

   if (findEdgeOut (idxFrom, idxTo) >= 0)  {
      return 1;
   }

   ssize_t pathLen = 0;
   SxSet<ssize_t> visited;
   SxList<SxPair<ssize_t,ssize_t> > distMap;

   visited << idxFrom;
   nextLevel (idxFrom, pathLen, maxPathLen, selection, &distMap, &visited);

   while (distMap.getSize () > 0)  {
      ssize_t ptr = distMap.first().key;
      pathLen = distMap.first().value;
      distMap.removeFirst ();

      ssize_t idx = findEdgeOut (ptr, idxTo);
      if (idx >= 0)  {
         return pathLen + 1;
      }
      nextLevel (ptr, pathLen, maxPathLen, selection, &distMap, &visited);
   }

   return 0;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxList<typename SxGraph<N,E,GS,ItPair>::ConstIterator>
SxGraph<N,E,GS,ItPair>::getPath (const N &from, const N &to,
                                 ssize_t maxPathLen,
                                 const Selection &selection) const
{
   SX_TRACE ();
   return getPath (begin(from), begin(to), maxPathLen, selection);
}


template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxList<typename SxGraph<N,E,GS,ItPair>::ConstIterator>
SxGraph<N,E,GS,ItPair>::getPath (
      const typename SxGraph<N,E,GS,ItPair>::ConstIterator &from,
      const typename SxGraph<N,E,GS,ItPair>::ConstIterator &to,
      ssize_t                                              maxPathLen,
      const Selection                                      &selection) const
{
   SX_TRACE ();
   SxList<ConstIterator> res;

   ssize_t idxFrom = from.getIdx();
   if (idxFrom < 0 || (selection.getPtr () && !selection->contains (idxFrom)))
      return res;
   ssize_t idxTo = to.getIdx();
   if (idxTo < 0 || (selection.getPtr () && !selection->contains (idxTo)))
      return res;

   if (findEdgeOut (idxFrom, idxTo) >= 0)  {
      res << ConstIterator (sx::Forward, this, idxFrom, Selection());
      res << ConstIterator (sx::Forward, this, idxTo, Selection());
   }  else  {
      ssize_t pathLen = 0;
      SxSet<ssize_t> visited;
      SxList<SxPair<ssize_t,ssize_t> > distMap;

      visited << idxFrom;
      nextLevel (idxFrom, pathLen, maxPathLen, selection, &distMap, &visited);

      while (distMap.getSize () > 0)  {
         ssize_t ptr = distMap.first().key;
         pathLen = distMap.first().value;
         distMap.removeFirst ();

         ssize_t idx = findEdgeOut (ptr, idxTo);
         if (idx >= 0)  {
            // --- reverse iteration to get the elements in path
            res << ConstIterator (sx::Forward, this, idxTo, Selection());
            visited.removeElement (idxTo);
            while (ptr != idxFrom)  {
               res.prepend (ConstIterator (sx::Forward, this, ptr,
                            Selection()));
               visited.removeElement (ptr);
               // --- parent

               const SxArray<SxPair<ssize_t,ssize_t> > &edgesIn = nodes->get(ptr).in ();
               ssize_t n = edgesIn.size;
               bool found = false;
               for (ssize_t i=0; i < n; ++i)  {
                  if (visited.contains (edgesIn(i).value))  {
                     ptr = edgesIn(i).value;
                     found = true;
                     break;
                  }
               }

               SX_CHECK (found); // inf loop
            }
            res.prepend (ConstIterator (sx::Forward, this, ptr, Selection()));
            break;
         }
         nextLevel (ptr, pathLen, maxPathLen, selection, &distMap, &visited);
      }
   }
   return res;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::nextLevel (ssize_t ptr, ssize_t pathLen,
                                        ssize_t maxPathLen,
                                        const Selection &selection,
                                        SxList<SxPair<ssize_t,ssize_t> > *distMap,
                                        SxSet<ssize_t> *visited) const
{
   SX_TRACE ();
   if (maxPathLen == 0 || pathLen + 1 < maxPathLen)  {
      SX_CHECK (ptr >= 0);
      SX_CHECK (distMap);
      SX_CHECK (visited);

      // --- search in children
      const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get(ptr).out ();
      if (selection.getPtr ())  {
         for (ssize_t i = 0; i < edgesOut.size; ++i)  {
            ssize_t key = edgesOut.elements[i].value;
            if (selection->contains (key)
                && !visited->contains (key))  {

               distMap->append (SxPair<ssize_t,ssize_t>(key,
                                                        pathLen + 1));
               visited->append (key);
            }
         }
      }  else  {
         for (ssize_t i = 0; i < edgesOut.size; ++i)  {
            ssize_t key = edgesOut.elements[i].value;
            if (!visited->contains (key))  {
               distMap->append (SxPair<ssize_t,ssize_t>(key,
                                                        pathLen + 1));
               visited->append (key);
            }
         }
      }

   }
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxList<N> SxGraph<N,E,GS,ItPair>::topsort (const Selection &selection_) const
{
   SX_TRACE ();
   SxList<N> res;
   SxList<ssize_t> s;
   ssize_t ptr = 0;

   // --- number of incoming edges for each vertex and start nodes

   ssize_t len = nodes->getSize ();
   SxArray<ssize_t> edges_(len);
   edges_.set (0);
   if (selection_.getPtr ())  {
      typename SxList<ssize_t>::ConstIterator it;
      for (it = selection_->begin(); it != selection_->end(); ++it)  {
         ssize_t nSelected = 0;
         const SxArray<SxPair<ssize_t,ssize_t> > &edgesIn = nodes->get(*it).in ();
         for (ssize_t i = 0; i < edgesIn.size; ++i)  {
            if (selection_->contains (edgesIn(i).value))  {
               nSelected++;
            }
         }
         edges_(*it) = nSelected;
         if (nSelected < 1)  s << (*it);
      }
   }  else  {
      ptr = nodes->first ();
      while (ptr >= 0)  {
         edges_(ptr) = nodes->get (ptr).in ().size;
         if (nodes->get (ptr).in ().size < 1)  s << ptr;
         ptr = nodes->next (ptr);
      }
   }

   // --- remove incoming edges and collect new start nodes
   while (s.getSize () > 0)  {
      ptr = s.first ();
      s.removeFirst ();
      res << nodes->get (ptr).element;

      const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get(ptr).out ();
      for (ssize_t i = 0; i < edgesOut.size; ++i)  {
         ssize_t idxTo = edgesOut(i).value;
         if (edges_(idxTo) > 0)  {
            edges_(idxTo)--;
            if (edges_(idxTo) < 1)  s << idxTo;
         }
      }
   }

   // --- check cycles
   for (ssize_t i = 0; i < len; ++i)  {
      if (edges_(i) > 0)  {
         SX_THROW ("graph has a cycle");
      }
   }
   return res;
}

// Enumerating Circuits and Loops in Graphs with Self-Arcs and Multiple-Arcs
// K.A. Hawick and H.A. James
// Proc. 2008 Int. Conf. on Foundations of Computer Science (FCS’08), p 14-20
// 2008
template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxList<SxList<ssize_t> > SxGraph<N,E,GS,ItPair>::getCyclesIdx (
   const Selection &selection_) const
{
   SX_TRACE ();
   SxList<SxList<ssize_t> > circuits;
   SxArray<bool> visited;
   SxArray<bool> blocked;
   SxArray<SxArray<ssize_t> > B; // recursive block/unblock
   SxList<ssize_t> stack;            // curent circuit

   ssize_t len = nodes->getSize ();
   blocked.resize (len);
   B.resize (len);
   visited.resize (len);

   // --- test circuit in each node
   if (selection_.getPtr ())  {
      if (len > 0)  visited.set (true);
      typename SxList<ssize_t>::ConstIterator it;
      for (it = selection_->begin(); it != selection_->end(); ++it)  {
         visited(*it) = false;
      }
      for (it = selection_->begin(); it != selection_->end(); ++it)  {
         blocked.set (false);
         for (ssize_t i=0; i < len; i++)  {
            if (B(i).getSize () > 0) B(i).removeAll ();
         }
         circuit (*it, *it, stack, visited, blocked, B, &circuits);
         visited(*it) = true;
      }
   }  else  {
      if (len > 0)  visited.set (false);
      ssize_t ptr = nodes->first ();
      while (ptr >= 0)  {
         blocked.set (false);
         for (ssize_t i=0; i < len; i++)  {
            if (B(i).getSize () > 0) B(i).removeAll ();
         }
         circuit (ptr, ptr, stack, visited, blocked, B, &circuits);
         visited(ptr) = true;
         ptr = nodes->next (ptr);
      }
   }
   return circuits;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxArray<SxArray<N> > SxGraph<N,E,GS,ItPair>::getCycles (
   const Selection &selection_) const
{
   SX_TRACE ();
   SxList<SxList<ssize_t> > circuits = getCyclesIdx (selection_);

   // --- convert graph indices to elements in found circuits
   SxArray<SxArray<N> > res;
   res.resize (circuits.getSize ());
   SxList<SxList<ssize_t> >::ConstIterator it;
   ssize_t i = 0;
   for (it = circuits.begin(); it != circuits.end(); ++it)  {
      res(i).resize (it->getSize ());
      SxList<ssize_t>::ConstIterator it2;
      ssize_t j = 0;
      for (it2 = it->begin(); it2 != it->end(); ++it2)  {
         res(i)(j) = nodes->get (*it2).element;
         j++;
      }
      i++;
   }
   return res;
}

// getCycles
template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
bool SxGraph<N,E,GS,ItPair>::circuit (ssize_t idx, ssize_t start,
                                      SxList<ssize_t> &stack,
                                      SxArray<bool> &visited,
                                      SxArray<bool> &blocked,
                                      SxArray<SxArray<ssize_t> > &B,
                                      SxList<SxList<ssize_t> > *circuits) const
{
   SX_TRACE ();
   SX_CHECK (circuits);

   bool res = false;
   stack.append (idx);
   blocked(idx) = true;


   const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get (idx).out ();
   for (ssize_t i = 0; i < edgesOut.size; ++i)  {
      ssize_t to = edgesOut(i).value;
      if (!visited(to))  {
         if (to == start)  {
            // --- new circuit
            // A -> A, A -> B -> A, A -> B -> C -> A, ...
            circuits->append (stack);
            res = true;
         }  else if (!blocked(to))  {
            if (circuit (to, start, stack, visited, blocked, B, circuits))  {
               res = true;
            }
         }
      }
   }

   if (res)  {
      unblock (idx, blocked, B);
   }  else  {
      for (ssize_t i = 0; i < edgesOut.size; ++i)  {
         ssize_t to = edgesOut(i).value;
         if (!visited(to) && !B(to).contains (idx))  {
            B(to).append (idx);
         }
      }
   }

   stack.removeLast ();
   return res;
}

// getCycles
template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::unblock (ssize_t idx,
                                      SxArray<bool> &blocked,
                                      SxArray<SxArray<ssize_t> > &B) const
{
   SX_TRACE ();
   blocked(idx) = false;

   for (ssize_t i = 0; i < B(idx).size; ++i)  {
      ssize_t e = B(idx)(i);
      B(idx).remove (i);
      //B(idx).remove (e);
      i--;
      if (blocked(e))  {
         unblock (e, blocked, B);
      }
   }
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::EdgeEvent &
SxGraph<N,E,GS,ItPair>::getSignal (const Iterator &it, sx::GraphEvent e)
{
   using namespace sx;
   ssize_t idx = it.getIdx ();
   switch (e)  {
      case BeforeEdgeInEvent  : return sigBeforeEdgeIn(idx);
      case BeforeEdgeOutEvent : return sigBeforeEdgeOut(idx);
      case AfterEdgeInEvent   : return sigAfterEdgeIn(idx);
      case AfterEdgeOutEvent  : return sigAfterEdgeOut(idx);
   }
   SX_EXIT;
   return sigBeforeEdgeOut (0); // dummy to make it compile
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::handleEdgeEvent (unsigned int e,
                                              ssize_t idxFrom,
                                              ssize_t idxTo)
{
   Selection sel;
   Iterator i1 (sx::Forward, this, idxFrom, sel);
   Iterator i2 (sx::Forward, this, idxTo,   sel);

   if (e & sx::BeforeEdgeInEvent && sigBeforeEdgeIn.containsKey (idxTo))
      sigBeforeEdgeIn(idxTo).send (i1, i2);
   if (e & sx::BeforeEdgeOutEvent && sigBeforeEdgeOut.containsKey (idxFrom))
      sigBeforeEdgeOut(idxFrom).send (i1, i2);
   if (e & sx::AfterEdgeInEvent && sigAfterEdgeIn.containsKey (idxTo))
      sigAfterEdgeIn(idxTo).send (i1, i2);
   if (e & sx::AfterEdgeOutEvent && sigAfterEdgeOut.containsKey (idxFrom))
      sigAfterEdgeOut(idxFrom).send (i1, i2);
}



template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::print () const
{

   SxUniqueList<ssize_t> visited;
   for (ssize_t i = 0; i < nodes->getSize (); ++i)  {
      const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = nodes->get (i).out ();
      for (ssize_t j = 0; j < edgesOut.getSize (); ++j)  {
         std::cout << nodes->get (i).element
                   << " --> "
                   << nodes->get (edgesOut.elements[j].value).element
                   << "\n";
      }
   }

}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void
SxGraph<N,E,GS,ItPair>::print (typename SxGraph<N,E,GS,ItPair>::Iterator it,
                               std::ostream &s,
                               const SxString &opt) const
{
   SxUniqueList<ssize_t> visited;
   typename SxGraph<N,E,GS,ItPair>::Iterator itEnd;
   for (; it != itEnd; ++it)  {
      visited << it.getIdx ();
   }
   print (visited, s, opt);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void
SxGraph<N,E,GS,ItPair>::print (typename SxGraph<N,E,GS,ItPair>::ConstIterator it,
                               std::ostream &s,
                               const SxString &opt) const
{
   SxUniqueList<ssize_t> visited;
   typename SxGraph<N,E,GS,ItPair>::ConstIterator itEnd;
   for (; it != itEnd; ++it)  {
      visited << it.getIdx ();
   }
   print (visited, s, opt);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
void SxGraph<N,E,GS,ItPair>::print (const SxUniqueList<ssize_t> &visited,
                                    std::ostream &s, const SxString &opt) const
{
   s << "digraph G {" << endl;
   if (opt != "")  {
      s << opt << endl;
   }

   typename SxList<ssize_t>::ConstIterator itIdx;
   for (itIdx = visited.begin(); itIdx != visited.end(); ++itIdx)  {
      const NT *node = &nodes->get (*itIdx);
      if (node->out ().size < 1 && node->in ().size < 1)  {
         s << "   \"" << node->element << "\"" << endl;
      }  else  {
         const SxArray<SxPair<ssize_t,ssize_t> > &edgesOut = node->out ();
         bool someEdgeVisited = false;
         for (ssize_t i = 0; i < edgesOut.size; ++i)  {
            if (visited.contains (edgesOut(i).value))  {
               s << "   \"" << node->element << "\" -> \""
                 << nodes->get (edgesOut(i).value).element
                 << "\"" << endl;
               someEdgeVisited = true;
            }
         }
         if (!someEdgeVisited && edgesOut.getSize () > 0)
            s << "   \"" << node->element << "\"" << endl;
      }
   }
   s << "}" << endl;

}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
ssize_t SxGraph<N,E,GS,ItPair>::findPos (const N &elem) const
{
   SX_TRACE ();
   return nodes->findPos (elem);
}


template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
size_t SxGraph<N,E,GS,ItPair>::getNBytes () const
{
   SX_TRACE ();
   size_t nBytes = sizeof(*this)
                 + nodes->getNBytes () + edges->getNBytes ();

   return nBytes;
}


template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
size_t getNBytes (const SxGraph<N,E,GS,ItPair> &in)
{
   return in.getNBytes ();
}
// ---------------------------------------------------------------------------

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::begin (ssize_t nodeIdx, sx::Direction dir_)
{
   SX_TRACE ();
   if (!nodes->contains(nodeIdx))
      nodeIdx = -1;
   return typename SxGraph<N,E,GS,ItPair>::Iterator (dir_, this, nodeIdx,
                                                     Selection());
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::begin (ssize_t nodeIdx, sx::Direction dir_) const
{
   SX_TRACE ();
   if (!nodes->contains(nodeIdx))
      nodeIdx = -1;
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator (dir_, this,
                                                          nodeIdx,
                                                          Selection());
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::begin (const Selection &selection_)
{
   SX_TRACE ();
   ssize_t idx = -1;
   if (selection_.getPtr ())  {
      if (selection_->getSize () > 0) idx = selection_->first ();
   }  else  {
      idx = nodes->first ();
   }

   return typename SxGraph<N,E,GS,ItPair>::Iterator (sx::Undefined, this,
                                                     idx, selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::begin (const Selection &selection_) const
{
   SX_TRACE ();
   ssize_t idx = -1;
   if (selection_.getPtr ())  {
      if (selection_->getSize () > 0) idx = selection_->first ();
   }  else  {
      idx = nodes->first ();
   }
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator (sx::Undefined,
                                                          this, idx,
                                                          selection_);

}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::begin (sx::Direction dir_, const N &elem_,
                               const Selection &selection_)
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::Iterator (dir_, this, idx,
                                                     selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::begin (sx::Direction dir_, const N &elem_,
                               const Selection &selection_) const
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator (dir_, this, idx,
                                                          selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::begin (const N &elem_, const Selection &selection_)
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::Iterator (sx::Forward, this, idx,
                                                     selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::begin (const N &elem_,
                               const Selection &selection_) const
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator(sx::Forward, this,
                                                         idx, selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::beginIn (const N &elem_, const Selection &selection_)
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::Iterator(sx::Backward, this, idx,
                                                    selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::beginBoth (const N &elem_,
                                   const Selection &selection_)
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::Iterator(sx::Both, this, idx,
                                                    selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::beginIn  (const N &elem_,
                                  const Selection &selection_) const
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator(sx::Backward, this,
                                                         idx, selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::beginBoth (const N &elem_,
                                   const Selection &selection_) const
{
   SX_TRACE ();
   ssize_t idx = findPos (elem_);
   if (selection_.getPtr () && !selection_->contains (idx)) idx = -1;
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator(sx::Both, this,
                                                         idx, selection_);
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::end ()
{
   return typename SxGraph<N,E,GS,ItPair>::Iterator ();
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::end () const
{
   return typename SxGraph<N,E,GS,ItPair>::ConstIterator ();
}
template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::Iterator
SxGraph<N,E,GS,ItPair>::getIterator (const typename SxGraph<N,E,GS,ItPair>::SelIdx &idx)
{
   SX_TRACE ();
   ssize_t ptr = idx;
   if (!nodes->contains (ptr))  ptr = -1;

   return typename SxGraph<N,E,GS,ItPair>::Iterator (sx::Forward, this,
                                                     ptr, Selection());
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::getIterator (const typename SxGraph<N,E,GS,ItPair>::SelIdx &idx) const
{
   SX_EXIT;
   SX_TRACE ();
   ssize_t ptr = idx;
   if (!nodes->contains (ptr))  ptr = -1;

   return typename SxGraph<N,E,GS,ItPair>::ConstIterator (sx::Forward, this,
                                                          ptr, Selection());
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
typename SxGraph<N,E,GS,ItPair>::ConstIterator
SxGraph<N,E,GS,ItPair>::getConstIterator (const typename SxGraph<N,E,GS,ItPair>::Iterator &it) const
{
   SX_EXIT;
   SX_TRACE ();
   return it;
}

template<class N,class E,
         template<class,bool> class GS,
         template<class,class,class> class ItPair>
SxPtr<SxGraph<N,E,GS,ItPair> >
SxGraph<N,E,GS,ItPair>::getContainer () const
{
   SX_EXIT;
   SX_TRACE ();
   return this->getThis ();
}

