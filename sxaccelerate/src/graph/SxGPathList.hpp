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
SxGPathList<N,E,GS>::SxGPathList ()
{
   // empty
}

template<class N,class E,template<class,bool> class GS>
SxGPathList<N,E,GS>::SxGPathList (const SxPtr<SxGraph<N,E,GS> > &gPtr_,
                                  const SelSet &sels_)
{
   SX_TRACE ();
   SX_CHECK (gPtr_.getPtr ());
   SX_CHECK (sels_.getPtr ());

   gPtr = gPtr_;
   sels = sels_;
}

template<class N,class E,template<class,bool> class GS>
SxGPathList<N,E,GS>::SxGPathList (const SxPtr<SxGraph<N,E,GS> > &gPtr_,
                                  const Selection &sel_)
{
   SX_TRACE ();
   SX_CHECK (gPtr_.getPtr ());
   SX_CHECK (sel_.getPtr ());

   gPtr = gPtr_;
   sels = SelSet::create ();
   sels->append (sel_);
}

template<class N,class E,template<class,bool> class GS>
SxGPathList<N,E,GS>::~SxGPathList ()
{
   // empty
}

template<class N,class E,template<class,bool> class GS>
ssize_t SxGPathList<N,E,GS>::getSize () const
{
   SX_TRACE ();
   return sels->getSize ();
}

template<class N,class E,template<class,bool> class GS>
ssize_t SxGPathList<N,E,GS>::getPathSize () const
{
   SX_TRACE ();
   if (sels->getSize () > 0)  return sels->first()->getSize ();
   else                       return 0;
}

template<class N,class E,template<class,bool> class GS>
SxPtr<SxGraph<N,E,GS> > SxGPathList<N,E,GS>::getGraphPtr () const
{
   SX_TRACE ();
   return gPtr;
}

template<class N,class E,template<class,bool> class GS>
typename SxGPathList<N,E,GS>::SelSet SxGPathList<N,E,GS>::getSelSet () const
{
   SX_TRACE ();
   return sels;
}

template<class N,class E,template<class,bool> class GS>
SxGPath<N,E,GS>
SxGPathList<N,E,GS>::operator() (ssize_t selIdx)
{
   SX_TRACE ();
   SX_CHECK (selIdx < sels->getSize (), selIdx, sels->getSize ());
   return SxGPath<N,E,GS>(gPtr, (*sels)(selIdx));
}

template<class N,class E,template<class,bool> class GS>
typename SxGPathList<N,E,GS>::Iterator
SxGPathList<N,E,GS>::begin ()
{
   SX_TRACE ();
   return typename SxGPathList<N,E,GS>::Iterator (this, sels->begin ());
}

template<class N,class E,template<class,bool> class GS>
typename SxGPathList<N,E,GS>::ConstIterator
SxGPathList<N,E,GS>::begin () const
{
   return typename SxGPathList<N,E,GS>::ConstIterator (this, sels->begin ());
}

template<class N,class E,template<class,bool> class GS>
typename SxGPathList<N,E,GS>::Iterator
SxGPathList<N,E,GS>::begin (ssize_t idx)
{
   SX_TRACE ();
   SX_CHECK (idx < sels->getSize (), idx, sels->getSize ());
   return typename SxGPathList<N,E,GS>::Iterator (this, sels->begin (idx));
}

template<class N,class E,template<class,bool> class GS>
typename SxGPathList<N,E,GS>::ConstIterator
SxGPathList<N,E,GS>::begin (ssize_t idx) const
{
   SX_TRACE ();
   SX_CHECK (idx < sels->getSize (), idx, sels->getSize ());
   return typename SxGPathList<N,E,GS>::ConstIterator (this, sels->begin (idx));
}

template<class N,class E,template<class,bool> class GS>
typename SxGPathList<N,E,GS>::Iterator
SxGPathList<N,E,GS>::end ()
{
   SX_TRACE ();
   return typename SxGPathList<N,E,GS>::Iterator (this, sels->end ());
}

template<class N,class E,template<class,bool> class GS>
typename SxGPathList<N,E,GS>::ConstIterator
SxGPathList<N,E,GS>::end () const
{
   SX_TRACE ();
   return typename SxGPathList<N,E,GS>::ConstIterator (this, sels->end ());
}

template<class N,class E,template<class,bool> class GS>
template<class Fn>
void SxGPathList<N,E,GS>::foreach (Fn fn)
{
   SX_TRACE ();
   sx::foreach (begin (), end (), fn);
}

template<class N,class E,template<class,bool> class GS>
template<class Fn>
void SxGPathList<N,E,GS>::foreach (Fn fn) const
{
   SX_TRACE ();
   sx::foreach (begin (), end (), fn);
}

