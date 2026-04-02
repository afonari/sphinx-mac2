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

#ifndef _SX_GRAPH_STORAGE_H_
#define _SX_GRAPH_STORAGE_H_

#define SX_GRAPH_STORAGE(TT,VV)                                               \
   public:                                                                    \
      ssize_t add (const TT &elem)                                            \
      {                                                                       \
         return this->addElem (elem);                                         \
      }                                                                       \
      ssize_t add (TT &&elem)                                                 \
      {                                                                       \
         return this->addElem (std::move(elem));                              \
      }                                                                       \
      void remove (ssize_t idx)                                               \
      {                                                                       \
         this->removeElem (idx);                                              \
      }                                                                       \
      VV &get (ssize_t idx)                                                   \
      {                                                                       \
         return this->getElem (idx);                                          \
      }                                                                       \
      const VV &get (ssize_t idx) const                                       \
      {                                                                       \
         return this->getElem (idx);                                          \
      }                                                                       \
      ssize_t findPos (const TT &elem)                                        \
      {                                                                       \
         return this->getIdx (elem);                                          \
      }                                                                       \
      ssize_t first () const                                                  \
      {                                                                       \
         return this->getFirstIdx ();                                         \
      }                                                                       \
      ssize_t next (ssize_t currIdx) const                                    \
      {                                                                       \
         return this->getNextIdx (currIdx);                                   \
      }                                                                       \
      ssize_t prev (ssize_t currIdx) const                                    \
      {                                                                       \
         return this->getPrevIdx (currIdx);                                   \
      }                                                                       \
      bool contains (ssize_t idx)                                             \
      {                                                                       \
         return this->containsElem (idx);                                     \
      }                                                                       \
      ssize_t getSize () const                                                \
      {                                                                       \
         return this->getNElems ();                                           \
      }                                                                       \
      size_t getNBytes () const                                               \
      {                                                                       \
         return this->getSizeBytes ();                                        \
      }                                                                       \
      void removeAll () {                                                     \
         this->removeElems ();                                                \
      }

#endif /*_SX_GRAPH_STORAGE_H_*/
