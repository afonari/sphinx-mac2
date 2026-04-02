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

#include <SxGQExprList.h>

SxGQExprList::SxGQExprList () : SxGQExprBase(SxGQExprBase::ExprType::Sibling),
                                sibType(OrderedDirect)
{
   // empty
}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                            SxGQExprList::SiblingType sType)
   : SxGQExprBase(SxGQExprBase::ExprType::Sibling),
     sibType(sType)
{
   SX_TRACE ();
   exprList << e1;
}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                            const SxPtr<SxGQExprBase> &e2,
                            SxGQExprList::SiblingType sType)
   : SxGQExprBase(SxGQExprBase::ExprType::Sibling),
     sibType(sType)
{
   SX_TRACE ();
   exprList << e1 << e2;
}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                            const SxPtr<SxGQExprList> &e2,
                            SxGQExprList::SiblingType sType)
   : SxGQExprBase(SxGQExprBase::ExprType::Sibling),
     sibType(sType)
{
   SX_TRACE ();
   SX_CHECK (sibType == e2->sibType);
   exprList << e1 << e2->exprList;
}

SxGQExprList::SxGQExprList (const SxPtr<SxGQExprList> &e1,
                            const SxPtr<SxGQExprBase> &e2,
                            SxGQExprList::SiblingType sType)
   : SxGQExprBase(SxGQExprBase::ExprType::Sibling),
     sibType(sType)
{
   SX_TRACE ();
   SX_CHECK (sibType == e1->sibType);
   exprList << e1->exprList << e2;
}

SxGQExprList::~SxGQExprList ()
{
   // empty
}

void SxGQExprList::append (const SxGQExprList &lst)
{
   SX_TRACE ();
   exprList << lst.exprList;
}

void SxGQExprList::append (const SxPtr<SxGQExprBase> &e)
{
   SX_TRACE ();
   SX_CHECK (e.getPtr ());
   exprList << e;
}

void SxGQExprList::append (const SxPtr<SxGQExprList> &lst)
{
   SX_TRACE ();
   SX_CHECK (lst.getPtr ());
   append (*lst);
}

SxGQExprList::SiblingType SxGQExprList::getSibType () const
{
   SX_TRACE ();
   return sibType;
}

const SxList<SxPtr<SxGQExprBase> > &SxGQExprList::getExprList () const
{
   SX_TRACE ();
   return exprList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExprList::firsts () const
{
   SX_TRACE ();
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   for (auto listIt = exprList.begin ();listIt != exprList.end (); ++listIt)
   {
      eList->append (*listIt);
   }
   return eList;
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQExprList::lasts () const
{
   SX_TRACE ();
   SxPtr<SxList<SxPtr<SxGQExprBase> > > eList =
      SxPtr<SxList<SxPtr<SxGQExprBase> > >::create ();

   for (auto listIt = exprList.begin ();listIt != exprList.end (); ++listIt)
   {
      if ((*listIt)->isOp ())  eList->append ((*listIt)->last ());
      else                     eList->append (*listIt);
   }
   return eList;
}

void SxGQExprList::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_TRACE ();
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   SX_CHECK (exprList.getSize () == lst->getSize());
   auto lastsIt = lst->begin ();
   for (auto listIt = exprList.begin();listIt != exprList.end(); ++listIt)
   {
      if ((*listIt)->isOp ())  (*listIt)->setLast (*lastsIt);
      else                     *listIt = *lastsIt;
      ++lastsIt;
   }
}

SxPtr<SxGQExprBase> SxGQExprList::first () const
{
   SX_TRACE ();
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   if (exprList.first ()->isOp())  return exprList.first ()->first ();
   else                            return exprList.first ();
}

SxPtr<SxGQExprBase> SxGQExprList::last () const
{
   SX_TRACE ();
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   if (exprList.last ()->isOp ())  return exprList.last ()->last ();
   else                            return exprList.last ();
}

void SxGQExprList::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_TRACE ();
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   if (exprList.last ()->isOp ())  exprList.last ()->setLast (p);
   else                            exprList.last () = p;
}

bool SxGQExprList::isOp () const
{
   SX_TRACE ();
   return true;
}

SxGQExprBase::ExprType SxGQExprList::getRightOp () const
{
   SX_TRACE ();
   SX_CHECK (exprList.getSize () > 0, exprList.getSize ());
   return SxGQExprBase::ExprType::Sibling;
}

size_t SxGQExprList::getHash () const
{
   SX_EXIT;
   return 0;
}

void SxGQExprList::makeGraph (SxGraph<SxGQPattern> *gPtr) const
{
   SX_TRACE ();
   for (auto listIt = exprList.begin (); listIt != exprList.end (); ++listIt)
   {
      if ((*listIt)->isOp ())  {
         (*listIt)->makeGraph (gPtr);
      }  else  {
         gPtr->createNode ((*listIt)->getGraphNode ());
      }
   }
}

SxGQPattern SxGQExprList::getGraphNode () const
{
   SX_EXIT;
}

std::ostream &SxGQExprList::print (std::ostream &s) const
{
   SX_TRACE ();
   for (auto listIt = exprList.begin (); listIt != exprList.end (); ++listIt)
   {
      s << "\n";
      s << *(*listIt);
      if (!(*listIt)->isOp())  s << ";\n";
   }
   return s;
}

