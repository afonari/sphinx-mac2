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

#include <SxGQAnd.h>

SxGQAnd::SxGQAnd () : SxGQExprBase(SxGQExprBase::ExprType::AND)
{
   // empty
}
SxGQAnd::SxGQAnd (const SxPtr<SxGQExprBase> &l, const SxPtr<SxGQExprBase> &r)
   : SxGQExprBase(SxGQExprBase::ExprType::AND), left (l), right (r)
{
   // empty
}
SxGQAnd::~SxGQAnd ()
{
   // empty
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQAnd::firsts () const
{
   SX_EXIT;
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   eList->append (first ());
   return eList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQAnd::lasts () const
{
   SX_TRACE ();
   if (right->isOp ())  {
      return right->lasts ();
   }  else  {
      SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();
      eList->append (right);
      return eList;
   }
}

SxPtr<SxGQExprBase> SxGQAnd::first () const
{
   SX_TRACE ();
   SX_CHECK (left.getPtr ());
   if (left->isOp())  return left->first ();
   else               return left;
}

SxPtr<SxGQExprBase> SxGQAnd::last () const
{
   SX_TRACE ();
   SX_CHECK (right.getPtr ());
   if (right->isOp ())  return right->last ();
   else                 return right;
}

void SxGQAnd::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_TRACE ();
   SX_CHECK (right.getPtr ());
   if (right->isOp ())  right->setLast (p);
   else                 right = p;
}

void SxGQAnd::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_TRACE ();
   SX_CHECK (right.getPtr ());
   SX_CHECK (right->isOp ());
   right->setLasts (lst);
}

bool SxGQAnd::isOp () const
{
   SX_TRACE ();
   return false;
}

SxGQExprBase::ExprType SxGQAnd::getRightOp () const
{
   SX_TRACE ();
   return right->getRightOp ();
}

size_t SxGQAnd::getHash () const
{
   SX_TRACE ();
   return ( left->getHash ()
          + SxHashFunction::hash (SxString("&&"))
          + right->getHash () );
}

void SxGQAnd::makeGraph (SxGraph<SxGQPattern> *gPtr) const
{
   SX_TRACE ();
   SxPtr<SxGQAnd> expr = SxPtr<SxGQAnd>::create (*this);
   SxGQPattern n = getGraphNode ();
   gPtr->createNode (n);
}

SxGQPattern SxGQAnd::getGraphNode () const
{
   SX_TRACE ();
   SX_CHECK (!left->isOp ());
   SxPtr<SxGQAnd> expr = SxPtr<SxGQAnd>::create (*this);
   SxGQPattern n((ssize_t)getHash ());
   n.setExpr (expr);
   return n;
}

std::ostream &SxGQAnd::print (std::ostream &s) const
{
   SX_TRACE ();
   s << "\"";
   s << *left;
   s << "&&";
   s << *right;
   s << "\"";
   return s;
}
