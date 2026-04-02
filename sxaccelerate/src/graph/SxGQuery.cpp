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

#include <SxGQuery.h>

SxGQuery::SxGQuery ()
{
   // empty
}

SxGQuery::SxGQuery (const SxGQuery &in)
{
   SX_TRACE ();
   this->rootExpr = in.rootExpr;
   makeGraph ();
}

SxGQuery::SxGQuery (const sx::N &n)
{
   SX_TRACE ();
   SxPtr<SxGQExprBase> e = n;
   this->rootExpr = e;
   makeGraph ();
}

SxGQuery::SxGQuery (const SxPtr<SxGQExprBase> &r)
{
   SX_TRACE ();
   SX_CHECK (r.getPtr ());
   rootExpr = r;
   makeGraph ();
}

SxGQuery::SxGQuery (const SxPtr<SxGQExprList> &r)
{
   SX_TRACE ();
   SX_CHECK (r.getPtr ());
   rootExpr = r;
   makeGraph ();
}

SxGQuery::~SxGQuery ()
{
   // empty
}

SxGQuery &SxGQuery::operator= (const SxGQuery &in)
{
   SX_TRACE ();
   if (this == &in) return *this;
   this->rootExpr = in.rootExpr;
   patternGraph = SxGraph<SxGQPattern>();
   makeGraph ();
   return *this;
}

// --- operator= allows to assign a SxGQExprBase object to SxGQuery
SxGQuery &SxGQuery::operator= (const SxPtr<SxGQExprBase> &r)
{
   SX_TRACE ();
   SX_CHECK (r.getPtr ());
   this->rootExpr = r;
   patternGraph = SxGraph<SxGQPattern>();
   makeGraph ();
   return *this;
}

SxGQuery &SxGQuery::operator= (const SxPtr<SxGQExprList> &r)
{
   SX_TRACE ();
   SX_CHECK (r.getPtr ());
   this->rootExpr = r;
   patternGraph = SxGraph<SxGQPattern>();
   makeGraph ();
   return *this;
}

void SxGQuery::makeGraph ()
{
   SX_TRACE ();
   if (rootExpr->isOp ())  {
      rootExpr->makeGraph (&patternGraph);
   }  else  {
      patternGraph.createNode (rootExpr->getGraphNode ());
   }
}


// -----------------------------------------------------------------------------

std::ostream &operator<< (std::ostream &s, const SxGQuery &q)
{
   s << "\n\ndigraph G {\n";
   s << *(q.rootExpr);
   s << "\n}\n\n";
   return s;
}

// -----------------------------------------------------------------------------

// --- returns SxGQExprBase representing a binary AND condition
SxPtr<SxGQExprBase>
operator&& (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQAnd> expr1 = SxPtr<SxGQAnd>::create (e1, e2);
   return expr1;
}

// --- returns SxGQExprBase representing a binary OR condition
SxPtr<SxGQExprBase>
operator|| (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQOr> expr1 = SxPtr<SxGQOr>::create (e1, e2);
   return expr1;
}

// -----------------------------------------------------------------------------
// Ordered Direct

// --- joins two child nodes represented by two SxGQExprBase objects
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedDirect);
   return ptr;
}

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedDirect);
   return ptr;
}

// --- joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------
// Unordered Indirect

// --- joins two child nodes represented by two SxGQExpr objects
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedIndirect);
   return ptr;
}

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedIndirect);
   return ptr;
}

// --- joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------
// Unordered Direct

// --- joins two child nodes represented by two SxGQExpr objects
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedDirect);
   return ptr;
}

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprBase> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::UnorderedDirect);
   return ptr;
}

// --- joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprList> &e1,
            const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------
// Ordered Indirect

// --- joins two child nodes represented by two SxGQExpr objects
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedIndirect);
   return ptr;
}

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   e1->append (e2);
   return e1;
}

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   SxPtr<SxGQExprList> ptr = SxPtr<SxGQExprList>::create (e1, e2,
                             SxGQExprList::SiblingType::OrderedIndirect);
   return ptr;
}

// --- joins two child nodes represented by two SxGQExprList objects
SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprList> &e1,
           const SxPtr<SxGQExprList> &e2)
{
   e1->append (e2);
   return e1;
}

// -----------------------------------------------------------------------------

// --- joins SxGQExprBase(parent) with SxGQExprList(children) as SxGQRelations
SxPtr<SxGQExprBase>
operator^ (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprList> &elist)
{
   if (e1->isOp ())  {
      if (e1->exprType == SxGQExprBase::ExprType::Relation)  {
         if (e1->getRightOp () == SxGQExprBase::ExprType::Relation)  {
            SxPtr<SxGQExprBase> p = e1->last ();
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (*elist);
            SxPtr<SxGQRelations> ptr
               = SxPtr<SxGQRelations>::create (p, child);
            e1->setLast (ptr);
            return e1;
         }  else if (e1->getRightOp () == SxGQExprBase::ExprType::Sibling)  {
            SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
            auto lastsIt = lst->begin ();

            SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts
               = SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();

            for (; lastsIt != lst->end (); ++lastsIt)  {
               // --- deep copy the sibling list, to avoid any problems due to shared obj
               SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (*elist);
               SxPtr<SxGQRelations> ptr
                  = SxPtr<SxGQRelations>::create (*lastsIt, child);
               newLasts->append (ptr);
            }

            e1->setLasts (newLasts);
            return e1;
         }

      }  else if (e1->exprType == SxGQExprBase::ExprType::Sibling)  {
         SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
         auto lastsIt = lst->begin ();

         SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();

         for (; lastsIt != lst->end (); ++lastsIt)  {
            // --- deep copy the sibling list, to avoid any problems due to shared obj
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (*elist);
            SxPtr<SxGQRelations> ptr
               = SxPtr<SxGQRelations>::create (*lastsIt, child);
            newLasts->append (ptr);
         }

         e1->setLasts (newLasts);
         return e1;
      }  else  {
         // --- currently no other operator has high enough
         //     precedence to appear as parent
         SX_EXIT;
      }
   }  else  {
      SxPtr<SxGQRelations> ptr
         = SxPtr<SxGQRelations>::create (e1, elist);
      return ptr;
   }
   SX_EXIT;
   return e1;
}

// --- joins SxGQExprBase(parent) with SxGQExprBase(children) as SxGQRelations
SxPtr<SxGQExprBase>
operator^ (const SxPtr<SxGQExprBase> &e1,
           const SxPtr<SxGQExprBase> &e2)
{
   if (e1->isOp ())  {
      if (e1->exprType == SxGQExprBase::ExprType::Relation)  {
         if (e1->getRightOp () == SxGQExprBase::ExprType::Relation)  {
            SxPtr<SxGQExprBase> p = e1->last ();
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
            SxPtr<SxGQRelations> ptr
               = SxPtr<SxGQRelations>::create (p, child);
            e1->setLast (ptr);
            return e1;
         }  else if (e1->getRightOp () == SxGQExprBase::ExprType::Sibling)  {
            SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
            auto lastsIt = lst->begin ();

            SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts
               = SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();
            for (; lastsIt != lst->end (); ++lastsIt)  {
               SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
               SxPtr<SxGQRelations> ptr
                  = SxPtr<SxGQRelations>::create (*lastsIt, child);
               newLasts->append (ptr);
            }
            e1->setLasts (newLasts);
            return e1;
         }

      }  else if (e1->exprType == SxGQExprBase::ExprType::Sibling)  {
         SxPtr<SxList<SxPtr<SxGQExprBase> > > lst = e1->lasts ();
         auto lastsIt = lst->begin ();

         SxPtr<SxList<SxPtr<SxGQExprBase> > > newLasts =
         SxPtr<SxList<SxPtr<SxGQExprBase> > >::create();
         for (; lastsIt != lst->end (); ++lastsIt)  {
            SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
            SxPtr<SxGQRelations> ptr
               = SxPtr<SxGQRelations>::create (*lastsIt, child);
            newLasts->append (ptr);
         }
         e1->setLasts (newLasts);
         return e1;
      }  else  {
         // --- currently no other operator has high enough
         //     precedence to appear as parent
         SX_EXIT;
      }
   }  else  {
      SxPtr<SxGQExprList> child = SxPtr<SxGQExprList>::create (e2);
      SxPtr<SxGQRelations> ptr = SxPtr<SxGQRelations>::create (e1, child);
      return ptr;
   }
   SX_EXIT;
   return e1;
}

