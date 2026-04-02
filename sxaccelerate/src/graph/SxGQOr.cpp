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

#include <SxGQOr.h>

SxGQOr::SxGQOr () : SxGQExprBase(SxGQExprBase::ExprType::OR)
{
   // empty
}
SxGQOr::SxGQOr (const SxPtr<SxGQExprBase> &l, const SxPtr<SxGQExprBase> &r)
   : SxGQExprBase(SxGQExprBase::ExprType::OR), left(l), right(r)
{
   // empty
}
SxGQOr::~SxGQOr ()
{
   // empty
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQOr::firsts () const
{
   SX_EXIT;
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   eList->append (first ());
   return eList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQOr::lasts () const
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

SxPtr<SxGQExprBase> SxGQOr::first () const
{
   SX_TRACE ();
   SX_CHECK (left.getPtr ());
   if (left->isOp())  return left->first ();
   else               return left;
}

SxPtr<SxGQExprBase> SxGQOr::last () const
{
   SX_TRACE ();
   SX_CHECK (right.getPtr ());
   if (right->isOp ())  return right->last ();
   else                 return right;
}

void SxGQOr::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_TRACE ();
   SX_CHECK (right.getPtr ());
   if (right->isOp ())  right->setLast (p);
   else                 right = p;
}

void SxGQOr::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_TRACE ();
   SX_CHECK (right.getPtr ());
   SX_CHECK (right->isOp ());
   right->setLasts (lst);
}

bool SxGQOr::isOp () const
{
   SX_TRACE ();
   return false;
}

SxGQExprBase::ExprType SxGQOr::getRightOp () const
{
   SX_TRACE ();
   return right->getRightOp ();
}


size_t SxGQOr::getHash () const
{
   SX_TRACE ();
   return ( left->getHash ()
          + SxHashFunction::hash (SxString("||"))
          + right->getHash () );
}

void SxGQOr::makeGraph (SxGraph<SxGQPattern> *gPtr) const
{
   SX_TRACE ();
   SxPtr<SxGQOr> expr = SxPtr<SxGQOr>::create (*this);
   SxGQPattern n = getGraphNode ();
   gPtr->createNode (n);
   gPtr->begin (n)->setExpr (expr);
}

SxGQPattern SxGQOr::getGraphNode () const
{
   SX_TRACE ();
   SX_CHECK (!left->isOp ());
   SxPtr<SxGQOr> expr = SxPtr<SxGQOr>::create (*this);
   SxGQPattern n((ssize_t)getHash ());
   n.setExpr (expr);
   return n;
}

std::ostream &SxGQOr::print (std::ostream &s) const
{
   SX_TRACE ();
   s << "\"";
   s << *left;
   s << "||";
   s << *right;
   s << "\"";
   return s;
}
