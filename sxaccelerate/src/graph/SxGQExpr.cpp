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

#include <SxGQExpr.h>

SxGQExpr::SxGQExpr ()
   : SxGQExprBase(SxGQExprBase::ExprType::None),
     captured(false),
     id(-1)
{
   // empty
}

SxGQExpr::SxGQExpr (const SxString &propName_,
                    const SxVariant &right_,
                    const SxGQExprBase::ExprType &type_)
   : SxGQExprBase(type_)
{
   SX_TRACE ();
   SX_CHECK (propName_!= "");
   propName = propName_;
   right    = right_;
   id       = (ssize_t)getHash ();
   captured = false;
}

SxGQExpr::SxGQExpr (const SxString &propName_,
                    const SxVariant &right_,
                    const SxString &capName_,
                    const SxGQExprBase::ExprType &type_,
                    bool isCap_)
   : SxGQExprBase(type_)
{
   SX_TRACE ();
   SX_CHECK (propName_!= "");
   propName     = propName_;
   captureName  = capName_;
   captured     = isCap_;
   right        = right_;
   id           = (ssize_t)getHash ();
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExpr::firsts () const
{
   SX_EXIT;
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   eList->append (first ());
   return eList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExpr::lasts () const
{
   SX_EXIT;
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   eList->append (last ());
   return eList;
}

SxPtr<SxGQExprBase> SxGQExpr::first() const
{
   SX_EXIT;
   SxPtr<SxGQExpr> e = SxPtr<SxGQExpr>::create (*this);
   return e;
}

SxPtr<SxGQExprBase> SxGQExpr::last() const
{
   SX_EXIT;
   SxPtr<SxGQExpr> e = SxPtr<SxGQExpr>::create (*this);
   return e;
}

void SxGQExpr::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_UNUSED (p);
   SX_EXIT;
}

void SxGQExpr::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_UNUSED (lst);
   SX_EXIT;
}

bool SxGQExpr::isOp () const
{
   SX_TRACE ();
   return false;
}

SxGQExprBase::ExprType SxGQExpr::getRightOp () const
{
   SX_TRACE ();
   return SxGQExprBase::ExprType::None;
}

bool SxGQExpr::isCaptured () const
{
   SX_TRACE ();
   return captured;
}

size_t SxGQExpr::getHash () const
{
   SX_TRACE ();
   size_t res = ( SxHashFunction::hash (propName)
                + SxHashFunction::hash (exprType) );
   if (right.isInitialized ())  res += SxHashFunction::hash (right.toString());
   return res;
}

void SxGQExpr::makeGraph (SxGraph<SxGQPattern> *gPtr) const
{
   SX_TRACE ();
   SxPtr<SxGQExpr> expr = SxPtr<SxGQExpr>::create (*this);
   SxGQPattern n = getGraphNode ();
   gPtr->createNode (n);
   gPtr->begin (n)->setExpr (expr);
}

SxGQPattern SxGQExpr::getGraphNode () const
{
   SX_TRACE ();
   SxPtr<SxGQExpr> expr = SxPtr<SxGQExpr>::create (*this);
   SxGQPattern n((ssize_t)getHash ());
   n.setExpr (expr);
   return n;
}

std::ostream &SxGQExpr::print (std::ostream &s) const
{
   SX_TRACE ();
   s << "\'N(" << propName << ")";
   if (exprType == SxGQExprBase::ExprType::NotEqual)  s << "!=";
   else                                               s << "==";
   if (exprType == SxGQExprBase::ExprType::Any)  s << "?";
   else                                          s << right.toString ();
   s << "\'";
   return s;
}
