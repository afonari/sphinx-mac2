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

template<class N,class E,template<class,bool> class GS>
SxGPath<N,E,GS>::SxGPath ()
{
   // empty
}

template<class N,class E,template<class,bool> class GS>
SxGPath<N,E,GS>::SxGPath (const SxPtr<SxGraph<N,E,GS> > &gPtr_,
                          const SxGPath::Selection &sel_)
{
   SX_TRACE ();
   SX_CHECK (gPtr_.getPtr ());
   SX_CHECK (sel_.getPtr ());

   gPtr = gPtr_;
   sel  = sel_;
}

template<class N,class E,template<class,bool> class GS>
SxGPath<N,E,GS>::~SxGPath ()
{
   // empty
}

template<class N,class E,template<class,bool> class GS>
ssize_t SxGPath<N,E,GS>::getSize () const
{
   SX_TRACE ();
   return sel->getSize ();
}

template<class N,class E,template<class,bool> class GS>
SxPtr<SxGraph<N,E,GS> > SxGPath<N,E,GS>::getGraphPtr () const
{
   SX_TRACE ();
   return gPtr;
}

template<class N,class E,template<class,bool> class GS>
typename SxGPath<N,E,GS>::Selection SxGPath<N,E,GS>::getSelection () const
{
   SX_TRACE ();
   return sel;
}

template<class N,class E,template<class,bool> class GS>
N &SxGPath<N,E,GS>::operator() (ssize_t idx)
{
   SX_TRACE ();
   ssize_t nodeIdx = -1;
   SX_CHECK (idx < sel->getSize (), idx, sel->getSize ());
   nodeIdx = (*sel)( idx );
   auto it = gPtr->begin (nodeIdx);
   SX_CHECK (it.isValid ());
   return *it;
}

template<class N,class E,template<class,bool> class GS>
template<class Fn>
void SxGPath<N,E,GS>::foreach (Fn fn)
{
   SX_TRACE ();
   SX_CHECK (gPtr.getPtr ());
   typename SxGraph<N,E,GS>::Iterator it = gPtr->begin (sel);
   for (; it != gPtr->end (); ++it)  {
      fn (it);
   }
}

template<class N,class E,template<class,bool> class GS>
template<class Fn>
void SxGPath<N,E,GS>::foreach (Fn fn) const
{
   SX_TRACE ();
   SX_CHECK (gPtr.getPtr ());
   typename SxGraph<N,E,GS>::ConstIterator it = gPtr->begin (sel);
   for (; it != gPtr->end (); ++it)  {
      fn (it);
   }
}
