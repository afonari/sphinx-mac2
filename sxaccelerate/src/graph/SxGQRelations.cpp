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
#include <SxGQRelations.h>

SxGQRelations::SxGQRelations ()
   : SxGQExprBase(SxGQExprBase::ExprType::Relation)
{
   // empty
}

SxGQRelations::SxGQRelations (const SxPtr<SxGQExprBase> &frm_,
                              size_t min_, size_t max_)
   : SxGQExprBase(SxGQExprBase::ExprType::Relation),
     from(frm_), minDist(min_), maxDist(max_)
{
   // empty
}

SxGQRelations::SxGQRelations (const SxPtr<SxGQExprBase> &frm_,
                              const SxPtr<SxGQExprList> &to_,
                              size_t min_, size_t max_)
   : SxGQExprBase(SxGQExprBase::ExprType::Relation),
     from(frm_), to(to_), minDist(min_),
     maxDist(max_)
{
   // empty
}

SxGQRelations::~SxGQRelations ()
{
   // empty
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQRelations::firsts () const
{
   SX_TRACE ();
   SX_CHECK (from.getPtr ());
   return from->firsts ();
}

SxPtr<SxList<SxPtr<SxGQExprBase> > > SxGQRelations::lasts () const
{
   SX_TRACE ();
   SX_CHECK (to.getPtr ());
   if (to->isOp ())  {
      return to->lasts ();
   }  else  {
      SxPtr<SxList<SxPtr<SxGQExprBase> > > res =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();
      res->append (to);
      return res;
   }
}

SxPtr<SxGQExprBase> SxGQRelations::first () const
{
   SX_TRACE ();
   SX_CHECK (from.getPtr ());
   if (from->isOp ())  return from->first ();
   else                return from;
}

SxPtr<SxGQExprBase> SxGQRelations::last () const
{
   SX_TRACE ();
   SX_CHECK (to.getPtr ());
   if (to->isOp ())  return to->last ();
   else              return to;
}

void SxGQRelations::setLast (const SxPtr<SxGQExprBase> &p)
{
   SX_TRACE ();
   SX_CHECK (to.getPtr ());
   if (to->isOp ())  to->setLast (p);
   else              to = p;
}

void SxGQRelations::setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst)
{
   SX_TRACE ();
   SX_CHECK (to.getPtr ());
   SX_CHECK (to->isOp ());
   to->setLasts (lst);
}

bool SxGQRelations::isOp () const
{
   SX_TRACE ();
   return true;
}

SxGQExprBase::ExprType SxGQRelations::getRightOp () const
{
   SX_TRACE ();
   return to->getRightOp ();
}


size_t SxGQRelations::getHash () const
{
   SX_EXIT;
   return 0;
}

void SxGQRelations::makeGraph (SxGraph<SxGQPattern> *gPtr) const
{
   SX_TRACE ();
   SxPtr<SxGQExprBase> last;
   if (from->isOp ())  {
      from->makeGraph (gPtr);
      last = from->last ();
   }  else  {
      last = from;
   }

   SxGQPattern parentNode  = last->getGraphNode ();

   gPtr->createNode (parentNode);

   // --- if node is already in graph, the new one will not be
   //     inserted, hence to make the changes to the node present
   //     in the graph use the returned iterator.

   gPtr->begin(parentNode)->setRelType (to->getSibType ());
   gPtr->begin(parentNode)->setExpr (last);

   for (auto listIt = to->getExprList().begin();
        listIt != to->getExprList().end(); ++listIt)
   {
      SxPtr<SxGQExprBase> child;
      if((*listIt)->isOp ())  child = (*listIt)->first ();
      else                    child = (*listIt);

      SxGQPattern childNode = child->getGraphNode ();

      gPtr->createNode (childNode);

      gPtr->begin (childNode)->setExpr (child);

      gPtr->createEdge (parentNode, childNode);

      if((*listIt)->isOp ())  (*listIt)->makeGraph (gPtr);
   }

}

SxGQPattern SxGQRelations::getGraphNode () const
{
   SX_EXIT;
}

std::ostream &SxGQRelations::print (std::ostream &s) const
{
   SX_TRACE ();
   if (from->isOp ())  s << *from;

   for (auto listIt = to->getExprList().begin();
        listIt != to->getExprList().end(); ++listIt)
   {
      s << "\n";

      if(from->isOp ())  s << *(from->last ());
      else               s << *from;

      s << " -> ";

      if((*listIt)->isOp ())  s << *((*listIt)->first ());
      else                    s << *(*listIt);

      s << ";\n";
      if((*listIt)->isOp ())  s << *(*listIt);
   }
   return s;
}

