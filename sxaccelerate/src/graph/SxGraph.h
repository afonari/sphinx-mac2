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

#ifndef _SX_GRAPH_H_
#define _SX_GRAPH_H_

#include <SxExportGraph.h>
#include <SxUniqueList.h>
#include <SxArrayList.h>
#include <SxHashedArray.h>
#include <SxIterator.h>
#include <SxContainer.h>
#include <SxGraphStorage.h>
#include <SxPair.h>
#include <SxSet.h>
#include <SxMap.h>
#include <SxTime.h>
#include <SxPtr.h>
#include <SxAlg.h>
#include <SxSignals.h>



namespace sx
{
   enum SX_EXPORT_GRAPH GraphEvent { BeforeEdgeInEvent  = 0x01,
                                     BeforeEdgeOutEvent = 0x02,
                                     AfterEdgeInEvent   = 0x04,
                                     AfterEdgeOutEvent  = 0x08 };
}

template<class T>
class SxGraphNode
{
   public:

      T element;

      ssize_t prev;
      ssize_t next;

      SxGraphNode () : prev(-1),next(-1),inE((ssize_t)0,(ssize_t)16),
                       outE((ssize_t)0,(ssize_t)16) { }

      inline SxArray<SxPair<ssize_t, ssize_t> > &in ()  { return outE; }
      inline SxArray<SxPair<ssize_t, ssize_t> > &out () { return inE;  }

      inline const SxArray<SxPair<ssize_t, ssize_t> > &in ()  const
      {
         return outE;
      }

      inline const SxArray<SxPair<ssize_t, ssize_t> > &out () const
      {
         return inE;
      }
   protected:
      SxArray<SxPair<ssize_t, ssize_t> > inE;
      SxArray<SxPair<ssize_t, ssize_t> > outE;
};

template<class T>
class SxGraphEdge
{
   public:

      T element;
      ssize_t prev;
      ssize_t next;
      SxGraphEdge (T val = T(), ssize_t in_ = -1, ssize_t out_ = -1,
                   ssize_t prev_ = -1, ssize_t next_ = -1)
                  : element (val), prev(prev_), next(next_),
                    inIdx (in_), outIdx (out_) { }

      ssize_t getOther (ssize_t v) {
         if (inIdx == v)        return outIdx;
         else if (outIdx == v)  return inIdx;
         else                   SX_EXIT; // invalid vertex
      }

      inline ssize_t &in  () { return inIdx;  }
      inline ssize_t &out () { return outIdx; }

      inline const ssize_t &in  () const { return inIdx;  }
      inline const ssize_t &out () const { return outIdx; }

   protected:
      ssize_t inIdx;
      ssize_t outIdx;
};


template<class T>
size_t getNBytes (const SxGraphNode<T> &in)
{
   size_t nBytes = sizeof (in)
                 + getNBytes (in.element) - sizeof (in.element)
                 + getNBytes (in.in ())   - sizeof (in.in ())
                 + getNBytes (in.out ())  - sizeof (in.out ());

   return nBytes;
}

template<class T>
size_t getNBytes (const SxGraphEdge<T> &in)
{
   size_t nBytes = sizeof (in)
                 + getNBytes (in.element) - sizeof (in.element);

   return nBytes;
}

// Provides a empty data type for edges.
// Contains auto incremented idx
// which results in unique hash value for each object.
class SX_EXPORT_GRAPH  SxBlankEdge
{
   public:
      SxBlankEdge () : idx(getIdx()) { }
      ~SxBlankEdge () { };
      operator ssize_t() const{ return idx; }
      bool operator== (const SxBlankEdge &in) { return idx == in.idx; }
      inline static ssize_t getIdx() { static ssize_t idx_(0); return idx_++; }

   protected:
      ssize_t idx;
};


template<class SValue,class SNode,class SContainer,class IT>
class SxGraphItState
{
   SX_ITERATOR_STATE
   // friend to State<const...>
   template<class CV,class CN,class CC,class CI> friend class SxGraphItState;

   public:

      using Selection = SxPtr<SxUniqueList<ssize_t> >;

      SxGraphItState ();
      SxGraphItState (sx::Direction dir_, SContainer *c, ssize_t idx_,
                      const Selection &sel);
      SxGraphItState (const SxGraphItState &in,
                      sx::ItCopyMode cMode_ = sx::CopyAll);
      SxGraphItState (SxGraphItState &&in, sx::ItCopyMode cMode_ = sx::CopyAll) noexcept;
      // non-const to const cast
      template<class CV,class CN,class CC,class CI>
         SxGraphItState (const SxGraphItState<CV,CN,CC,CI> &in_)
         : container(in_.container), node(in_.node) { copy(in_); }

      SxGraphItState &operator= (const SxGraphItState &);
      SxGraphItState &operator= (SxGraphItState &&);

      // --- SxGraph specifics
      Selection getUnvisited () const;


      IT in  (ssize_t) const;
      IT out (ssize_t) const;
      SValue &getIn (ssize_t);                 // deprecate, use in(ssize_t)
      SValue &getOut (ssize_t);                // deprecate, use out(ssize_t)
      ssize_t getSizeIn () const;
      ssize_t getSizeOut () const;
      ssize_t findEdgeOut (ssize_t idx);
      ssize_t findEdgeIn (ssize_t idx);

      inline ssize_t getDepth () const { return pathLen; }
      inline SNode *getNode () { return node; }
      inline ssize_t getIdx () const { return idx; }
      inline typename SContainer::SelIdx getSelIdx () const {
         SX_CHECK (container);
         return idx;
      }

      IT neighbors (ssize_t maxDepth=1) const;   // excludes parent
      IT hops (ssize_t maxDepth) const;   // includes parent
      IT reverse () const;


   protected:

      SContainer *container;   // graph
      SNode *node;
      Selection selection;
      ssize_t idx;
      ssize_t distMapIdx;
      ssize_t maxDepth;  // 0 if unlimited
      ssize_t pathLen;
      typename SxList<ssize_t>::ConstIterator selectionIt;
      //SxList<SxPair<ssize_t,ssize_t> > distMap;
      SxHashedArray<SxPair<ssize_t,ssize_t> > distMap;
      SxSet<ssize_t> visited;

      // --- selected neighbors
      mutable ssize_t selectedSizeIn;
      mutable ssize_t selectedSizeOut;
      mutable SxArray<ssize_t> selectedIn;
      mutable SxArray<ssize_t> selectedOut;

      void procNext ();
      void procPrev ();

      // --- SX_ITERATOR / SX_CONST_ITERATOR callbacks
      template<class CIT> void copy (const CIT &,
                                     sx::ItCopyMode cMode_ = sx::CopyAll);
      void move (SxGraphItState &&, sx::ItCopyMode cMode_ = sx::CopyAll);

      IT insertElem (ssize_t newPos, const SValue &elem);
      IT prependElem (const SValue &elem);
      IT appendElem (const SValue &elem);
      void next ();
      void prev ();
      bool valid () const { return (container != NULL && node != NULL); }
      SValue &getRef ();
      SValue *getPtr ();
      bool equal (const IT &in) const;

};

template<class T,class Node,class Container>
class SxGraphIterators
{
   public:
      using Selection = SxPtr<SxUniqueList<ssize_t> >;

      template<class V,class N,class C,class IT>
      using State = SxGraphItState<V,N,C,IT>;

      class Iterator : public State<T,Node,Container,Iterator>
      {
         SX_ITERATOR(T,Container,Iterator)
         public:

            using EdgeEvent = SxSignal<const Iterator&,
                                       const Iterator&,
                                       const char *>;


            Iterator ()
               : State<T,Node,Container,Iterator> () { }
            Iterator (sx::Direction dir_, Container *c, ssize_t idx_,
                      const Selection &sel)
               : State<T,Node,Container,Iterator> (dir_, c, idx_, sel)
            { }
            Iterator (const Iterator &in, sx::ItCopyMode cMode_ = sx::CopyAll)
               : State<T,Node,Container,Iterator> (in, cMode_)
                 { }
            Iterator (Iterator &&in, sx::ItCopyMode cMode_ = sx::CopyAll) noexcept
               : State<T,Node,Container,Iterator>
                 (std::move(in), cMode_) { }

            // conversion from State
            Iterator (const State<T,Node,Container,Iterator> &in_)
               : State<T,Node,Container,Iterator> (in_) { }
      };

      class ConstIterator
         : public State<const T,const Node,const Container,ConstIterator>
      {
         SX_CONST_ITERATOR(T,Container,ConstIterator)
         public:

            ConstIterator ()
               : State<const T,const Node,const Container,ConstIterator> () { }
            ConstIterator (sx::Direction dir_, const Container *c, ssize_t idx_,
                           const Selection &sel)
               : State<const T,const Node,const Container,ConstIterator> (
                    dir_, c, idx_, sel) { }
            ConstIterator (const ConstIterator &in,
                           sx::ItCopyMode cMode_ = sx::CopyAll)
               : State<const T,const Node,const Container,ConstIterator>
                 (in, cMode_) { }
            ConstIterator (ConstIterator &&in,
                           sx::ItCopyMode cMode_ = sx::CopyAll)
               : State<const T,const Node,const Container,ConstIterator>
                 (std::move(in), cMode_) { }

            // conversion from State
            ConstIterator (
               const State<const T,const Node,const Container,ConstIterator> &in_)
               : State<const T,const Node,const Container,ConstIterator> (in_)
                 { }

            // non-const to const cast
            ConstIterator (const Iterator &in_)
               : State<const T,const Node,const Container,ConstIterator>(in_)
                 { }

      };
};




template<class T,bool isNS>
class SxGraphStorage
{
   template<class N,class E,
            template<class,bool> class GS,
            template<class,class,class> class ItPair> friend class SxGraph;
   template<class CV,class CN,class CC,class CI> friend class SxGraphItState;


   protected:
      template<class N,bool Cond>
      class IsNodeStorage;

      template<class N>
      class IsNodeStorage<N,true>
      {
         public:
            typedef SxGraphNode<N> Elem;
      };

      template<class N>
      class IsNodeStorage<N,false>
      {
         public:
            typedef SxGraphEdge<N> Elem;
      };


      typedef typename IsNodeStorage<T,isNS>::Elem  Elem;

      SxHashedArray<Elem> elems;
      SxArray<ssize_t>    hashTable;

      ssize_t firstElement;
      ssize_t lastElement;
      ssize_t firstFree;

      ssize_t nElems;

      ssize_t hashSize;
      ssize_t hashUsed;
      ssize_t hashBound;


      // --- not to be called directly
      ssize_t     addElem (const T &);
      ssize_t     addElem (T &&);
      void        removeElem (ssize_t idx);

      bool    containsElem (ssize_t idx) const;
      ssize_t getNElems () const;
      size_t  getSizeBytes () const;
      void    removeElems ();

            Elem &getElem (ssize_t idx);
      const Elem &getElem (ssize_t idx) const;
      ssize_t     getIdx (const T &);
      ssize_t     getFirstIdx () const;
      ssize_t     getNextIdx  (ssize_t currIdx) const;
      ssize_t     getPrevIdx  (ssize_t currIdx) const;

      // --- hash table
      ssize_t hashFindPos (const T &) const;
      bool    hashContains (const T &, ssize_t *);
      void    hashResize (ssize_t newSize_);

   public:
      SX_GRAPH_STORAGE(T,Elem);

      SxGraphStorage ();
      SxGraphStorage (const SxGraphStorage<T,isNS> &);
      SxGraphStorage (SxGraphStorage<T,isNS> &&);
      SxGraphStorage<T,isNS> &operator= (const SxGraphStorage<T,isNS> &);
      SxGraphStorage<T,isNS> &operator= (SxGraphStorage<T,isNS> &&);
      ~SxGraphStorage ();

};


template<class N,class E=SxBlankEdge,
         template<class,bool> class GS=SxGraphStorage,
         template<class,class,class> class ItPair=SxGraphIterators>
class SxGraph : public SxThis<SxGraph<N,E,GS,ItPair> >
{
   public:
      template<class CV,class CN,class CC,class CI> friend class SxGraphItState;


      // --- Iterator
      using NT            = typename GS<N,true>::Elem;
      using ET            = typename GS<E,false>::Elem;
      using Container     = SxGraph<N,E,GS,ItPair>;
      using Iterator      = typename ItPair<N,NT,Container>::Iterator;
      using ConstIterator = typename ItPair<N,NT,Container>::ConstIterator;



      using Selection = SxPtr<SxUniqueList<ssize_t> >;


      // --- Event API
      using EdgeEvent = SxSignal<const Iterator&,
                                 const Iterator&,
                                 const char *>;

      // --- SxSelection:
      typedef N                     TElem;
      typedef ssize_t               SelIdx;
      typedef SxUniqueList<SelIdx>  SelStorage;
      typedef Container             SelContainer;

      SX_CONTAINER (Container)

      SxGraph ();
      SxGraph (const SxGraph<N,E,GS,ItPair> &);
      SxGraph (SxGraph<N,E,GS,ItPair> &&);
      virtual ~SxGraph ();

      SxGraph<N,E,GS,ItPair> &operator=(const SxGraph<N,E,GS,ItPair> &);
      SxGraph<N,E,GS,ItPair> &operator=(SxGraph<N,E,GS,ItPair> &&) noexcept;

      // -- Insertion
      virtual ConstIterator createNode (const N &);
      virtual ConstIterator createNode (N &&);
      virtual void createEdge (const N &from, const N &to, const E &elem=E());
      virtual void createEdge (ConstIterator &, ConstIterator &, const E &elem=E());

      // -- Deletion
      virtual void removeNode (ssize_t nodeIdx);
      virtual void removeEdge (const N &, const N &);
      void removeElement (const N &);
      void removeSelection (const Selection &selection);
      void removeAll ();

      // -- Search
      bool containsNode (const N &) const;
      ssize_t findPos (const N &) const;

      ssize_t findPath      (const N         &from,
                             const N         &to,
                             ssize_t         maxPathLen=0,
                             const Selection &selection=Selection()) const;

      SxList<ConstIterator>
         getPath (const N         &from,
                  const N         &to,
                  ssize_t         maxPathLen=0,
                  const Selection &selection=Selection()) const;

      SxList<ConstIterator>
         getPath (const ConstIterator &from,
                  const ConstIterator &to,
                  ssize_t             maxPathLen=0,
                  const Selection     &selection=Selection()) const;


      SxList<N> topsort (const Selection &selection=Selection()) const;

      SxArray<SxArray<N> > getCycles (const Selection &selection=Selection()) const;
      SxList<SxList<ssize_t> > getCyclesIdx(const Selection &selection=Selection()) const;


      // -- Size
      inline ssize_t getSize () const { return nodes->getSize (); }
      inline ssize_t getNEdges () const { return edges->getSize (); }
      size_t getNBytes () const;

      // Iteratorfunctions
      Iterator      begin (ssize_t nodeIdx, sx::Direction dir_=sx::Forward);
      ConstIterator begin (ssize_t nodeIdx, sx::Direction dir_=sx::Forward) const;

      Iterator      begin (const Selection &selection=Selection());
      ConstIterator begin (const Selection &selection=Selection()) const;
      Iterator      begin (sx::Direction dir_, const N      &from,
                           const Selection &selection=Selection());
      ConstIterator begin (sx::Direction dir_, const N      &from,
                           const Selection &selection=Selection()) const;
      Iterator      begin (const N &from,
                          const Selection &selection=Selection());
      ConstIterator begin (const N &from,
                           const Selection &selection=Selection()) const;
      Iterator beginIn   (const N &from,
                          const Selection &selection=Selection());
      Iterator beginBoth (const N &from,
                          const Selection &selection=Selection());
      Iterator end ();

      ConstIterator beginIn   (const N         &from,
                               const Selection &selection=Selection()) const;
      ConstIterator beginBoth (const N         &from,
                               const Selection &selection=Selection()) const;
      ConstIterator end () const;

      Iterator      fromLast ()  {SX_EXIT;}
      ConstIterator fromLast () const {SX_EXIT;}

      Iterator      toFirst ()  {SX_EXIT;}
      ConstIterator toFirst () const {SX_EXIT;}


      // --- SxSelection functions
      inline Iterator      getIterator (const SelIdx &);
      inline ConstIterator getIterator (const SelIdx &) const;
      inline ConstIterator getConstIterator (const Iterator &) const;
      SxPtr<Container> getContainer () const;

      // --- event handling
      EdgeEvent &getSignal (const Iterator &,sx::GraphEvent);

      void print (Iterator, std::ostream &,
                  const SxString &options="") const;
      void print (ConstIterator, std::ostream &,
                  const SxString &options="") const;

      void print () const;

      void print (const SxUniqueList<ssize_t> &visited,
                  std::ostream                &,
                  const SxString              &options="") const;

      // --- findPath
      void nextLevel (ssize_t                          ptr,
                      ssize_t                          pathLen,
                      ssize_t                          maxPathLen,
                      const Selection                  &selection,
                      SxList<SxPair<ssize_t,ssize_t> > *queue,
                      SxSet<ssize_t>                   *visited) const;

      // --- getCycles
      bool circuit (ssize_t                        idx,
                    ssize_t                        start,
                    SxList<ssize_t>                &stack,
                    SxArray<bool>                  &visited,
                    SxArray<bool>                  &blocked,
                    SxArray<SxArray<ssize_t> >     &B,
                    SxList<SxList<ssize_t> >       *circuits) const;
      void unblock (ssize_t, SxArray<bool> &, SxArray<SxArray<ssize_t> >&) const;


      bool hasOutEdge (ssize_t, ssize_t) const;
      ssize_t findEdgeOut (ssize_t, ssize_t) const;
      ssize_t findEdgeIn (ssize_t, ssize_t) const;
      void unlinkOut (ssize_t, ssize_t);
      void unlinkIn (ssize_t, ssize_t);



   protected:

      // --- events
      SxMap<ssize_t,EdgeEvent> sigBeforeEdgeIn;
      SxMap<ssize_t,EdgeEvent> sigBeforeEdgeOut;
      SxMap<ssize_t,EdgeEvent> sigAfterEdgeIn;
      SxMap<ssize_t,EdgeEvent> sigAfterEdgeOut;


      SxPtr<GS<N,true> > nodes;
      SxPtr<GS<E,false> > edges;

      void handleEdgeEvent (unsigned int e, ssize_t from, ssize_t to);


};


#include <SxGraph.hpp>

#endif /* _SX_GRAPH_H_ */
