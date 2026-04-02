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

#ifndef _SX_G_QUERY_H_
#define _SX_G_QUERY_H_

#include <SxGQExprBase.h>
#include <SxN.h>
#include <SxGQExpr.h>
#include <SxGQExprList.h>
#include <SxGQAnd.h>
#include <SxGQOr.h>
#include <SxGQRelations.h>


/** \brief Graph Query Class

    \b SxGQuery = SPHInX Graph Query Class

    SxGQuery is the top level class used to parse
    then store the given query. It is also used
    to find matching graph patterns.

 */
class SX_EXPORT_GRAPH SxGQuery
{

   public:
      typedef SxPtr<SxUniqueList<ssize_t> > Selection;
      typedef SxPtr<SxList<Selection> >  SelSet;

      SxGQuery ();
      SxGQuery (const SxGQuery &in);
      SxGQuery (const sx::N &n);
      SxGQuery (const SxPtr<SxGQExprBase> &r);
      SxGQuery (const SxPtr<SxGQExprList> &r);
     ~SxGQuery ();

      SxGQuery &operator= (const SxGQuery &in);
      SxGQuery &operator= (const SxPtr<SxGQExprBase> &r);
      SxGQuery &operator= (const SxPtr<SxGQExprList> &r);

      void makeGraph ();


      // --- returns the matched selections of nodes' Idx.
      template<class N,class E,template<class,bool> class GS>
      SelSet matchAll (const SxPtr<SxGraph<N,E,GS> > &g,
                       const typename SxGraph<N,E,GS>::Iterator &gIt);
      template<class N,class E,template<class,bool> class GS>
      SelSet matchAll (const SxPtr<SxGraph<N,E,GS> > &g);

      template<class N,class E,template<class,bool> class GS>
      Selection match (const SxPtr<SxGraph<N,E,GS> > &g,
                       const typename SxGraph<N,E,GS>::Iterator &gIt);
      template<class N,class E,template<class,bool> class GS>
      Selection  match (const SxPtr<SxGraph<N,E,GS> > &g);


      friend SX_EXPORT_GRAPH std::ostream &operator<< (std::ostream &s,
                                                       const SxGQuery &q);

   protected:
      SxPtr<SxGQExprBase> rootExpr;
      SxGraph<SxGQPattern> patternGraph;
      SelSet selections;

      inline void resizeSelections (const SelSet &sels,
                                    ssize_t setSize,
                                    ssize_t selSize) const;
      inline SelSet crossSelections (const SelSet &sels1,
                                     const SelSet &sels2) const;

      // ---------------------------------------------------------------------
      // match all functions

      /* \brief Find all matches for given node in incoming nodes

         findes all matches in incoming nodes of given node of
         data graph for the given node in pattern graph.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchAllInNode (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                           const SelSet &sels,
                           const Selection &visited) const;

      /* \brief Find all matches for all incoming nodes

         finds all matches of all nodes on the incoming
         edges of the pattern node, in the data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchAllInNodes (const SxGraph<SxGQPattern>::ConstIterator &it,
                            const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                            const SelSet &sels,
                            const Selection &visited) const;

      /* \brief Find all matches for given and connected nodes

         finds all matches for given pattern node and
         recursively goes through all connected nodes.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchAllRecursive (const SxGraph<SxGQPattern>::ConstIterator &it,
                              const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                              const SelSet &sels,
                              const Selection &visited) const;

      /* \brief Match all outgoing ordered-direct starting at childIdx

         finds all ordered-direct matches for nodes on outgoing
         edges of pattern node starting at given childIdx in
         data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchAllOutNodesODIdx (const SxGraph<SxGQPattern>::ConstIterator &it,
                                         const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                         ssize_t childIdx, const SelSet &sels,
                                         const Selection &visited) const;

      /* \brief Match all outgoing nodes ordered-direct

         finds all ordered-direct matches for nodes on outgoing
         edges of pattern node with given parent node in
         data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchAllOutNodesOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                      const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                      const SelSet &sels,
                                      const Selection &visited) const;

      /* \brief Match all outgoing nodes unordered-indirect

         finds all unordered-indirect matches for nodes on
         outgoing edges of pattern node with given parent
         node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchAllOutNodesUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                               const typename SxGraph<N,E,GS>::ConstIterator &pIt,
                               const SelSet &sels,
                               const Selection &visited) const;

      /* \brief Match all outgoing nodes ordered-indirect

         finds all ordered-indirect matches for nodes on
         outgoing edges of pattern node with given parent
         node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchAllOutNodesOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                      const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                      const SelSet &sels,
                                      const Selection &visited) const;

       /* \brief Match all outgoing nodes unordered-direct

          finds all unordered-direct matches for nodes on
          outgoing edges of pattern node with given parent
          node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchAllOutNodesUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                      const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                      const SelSet &sels,
                                      const Selection &visited) const;

      // ---------------------------------------------------------------------
      // single match functions

      /* \brief  Find one match of one node

         finds single match for given pattern node in
         the incoming nodes of data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchInNode (const SxGraph<SxGQPattern>::ConstIterator &it,
                        const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                        const Selection &sel,
                        const Selection &visited) const;

      /* \brief Find one match for all incoming nodes

         finds single match for all nodes on the incoming
         edges of the pattern node in incoming nodes of
         data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchInNodes (const SxGraph<SxGQPattern>::ConstIterator &it,
                         const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                         const Selection &sel,
                         const Selection &visited) const;

      /* \brief Find one match for given and connected nodes

         finds single match for given pattern node and
         recursively goes through all connected nodes.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchRecursive (const SxGraph<SxGQPattern>::ConstIterator &it,
                           const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                           const Selection &sel,
                           const Selection &visited) const;


      /* \brief Find an unordered-direct match starting at childIdx

         finds single unordered-direct match for nodes on
         outgoing edges of pattern node starting at a given
         childIdx in pattern and data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchOutNodesUDIdx (const SxGraph<SxGQPattern>::ConstIterator &it,
                               const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                               ssize_t childIdx, ssize_t exprChildIdx,
                               const Selection &sel,
                               const Selection &visited) const;

      /* \brief Find an unordered-direct match for all outgoing nodes

         finds single unordered-direct match for nodes on
         outgoing edges of pattern node with given parent
         node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchOutNodesUD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                   const Selection &sel,
                                   const Selection &visited) const;

      /* \brief Find an unordered-indirect match for outgoing nodes

         finds single unordered-indirect match for nodes on
         outgoing edges of pattern node with given parent
         node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchOutNodesUI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                   const Selection &sel,
                                   const Selection &visited) const;

      /* \brief Find an ordered-direct match starting at childIdx

         finds single ordered-direct match for nodes on
         outgoing edges of pattern node starting from
         given childIdx of parent node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      bool matchOutNodesODIdx (const SxGraph<SxGQPattern>::ConstIterator &pIt,
                               const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                               const Selection &sel,
                               ssize_t childIdx,
                               const Selection &visited) const;

      /* \brief Find an ordered-direct match for outgoing nodes

         finds single ordered-direct match for nodes on
         outgoing edges of pattern node with given parent
         node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchOutNodesOD (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                   const Selection &sel,
                                   const Selection &visited) const;

      /* \brief Find an ordered-indirect match for outgoing nodes

         finds single ordered-indirect match for nodes on
         outgoing edges of pattern node with given parent
         node in data graph.
      */
      template<class N,class E,template<class,bool> class GS>
      inline bool matchOutNodesOI (const SxGraph<SxGQPattern>::ConstIterator &it,
                                   const typename SxGraph<N,E,GS>::ConstIterator &gIt,
                                   const Selection &sel,
                                   const Selection &visited) const;
};

SX_EXPORT_GRAPH std::ostream &operator<< (std::ostream &s, const SxGQuery &q);

// -----------------------------------------------------------------------------

#include <SxGQuery.hpp>

// -----------------------------------------------------------------------------

// --- namespace to contain expression specific
//     match and matchAll functions.
namespace SxGQInternal  {

typedef SxPtr<SxUniqueList<ssize_t> > Selection;
typedef SxPtr<SxList<Selection> >  SelSet;

// --- finds one match for SxGQExpr
template<class N,class E,template<class,bool> class GS>
bool match (const SxPtr<SxGQExpr> &expr,
            const typename SxGraph<N,E,GS>::ConstIterator &it,
            const Selection &sel)
{
   SX_TRACE ();
   if (!it.isValid ()) return false;
   bool res = false;
   switch (expr->exprType)  {
      case SxGQExpr::ExprType::Any:
         if (it->hasProperty(expr->propName))  {
            res = true;
            sel->append (it.getIdx ());
         }
         return res;
      case SxGQExpr::ExprType::Equal:
         res =   (it->hasProperty (expr->propName)
               && it->getProperty(expr->propName) == expr->right);

         if (res == true) sel->append (it.getIdx ());
         if (res == false)
            SX_DBG_MSG ("SxGQExpr::Equal 'False' for Idx: "+
                        SxString(it.getIdx ())+
                        " == "+expr->right.toString ());
         return res;
      case SxGQExpr::ExprType::NotEqual:
         res =  !(it->hasProperty (expr->propName)
               && it->getProperty(expr->propName) == expr->right);
         if (res == true) sel->append (it.getIdx ());
         if (res == false)
            SX_DBG_MSG ("SxGQExpr::NotEqual 'False' for Idx: "+
                        SxString(it.getIdx ())+
                        " == "+expr->right.toString ());
         return res;
      default:
         return false;
   }
   return false;
}

// --- finds all matches for SxGQExpr
template<class N,class E,template<class,bool> class GS>
bool matchAll (const SxPtr<SxGQExpr> &expr,
               const typename SxGraph<N,E,GS>::ConstIterator &it,
               const SelSet &sels)
{
   SX_TRACE ();
   Selection sel = Selection::create ();
   if ( match<N,E,GS> (expr, it, sel) )  {

      SX_CHECK (sel->getSize () > 0, sel->getSize());

      // if empty, then add as first
      if (sels->getSize () == 0)  {
         sels->append (sel);
         return true;
      }

      // cross of n x 1 = n
      for (auto it1 = sels->begin ();it1 != sels->end (); ++it1)  {
         (*it1)->append (*sel);
      }
      return true;
   }
   return false;
}

// --- finds one match for SxGQAnd
template<class N,class E,template<class,bool> class GS>
bool match (const SxPtr<SxGQAnd> &expr,
            const typename SxGraph<N,E,GS>::ConstIterator &it,
            const Selection &sel)
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (expr->left.getPtr ());
   SX_CHECK (expr->right.getPtr ());

   bool res = match<N,E,GS> (expr->left, it, sel);
   if (res)  {
      Selection tmpSel = Selection::create ();
      tmpSel->append (*sel);
      res = res && match<N,E,GS> (expr->right, it, tmpSel);
   }
   return res;
}

// --- finds all matches for SxGQAnd
template<class N,class E,template<class,bool> class GS>
bool matchAll (const SxPtr<SxGQAnd> &expr,
               const typename SxGraph<N,E,GS>::ConstIterator &it,
               const SelSet &sels)
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (expr->left.getPtr ());
   SX_CHECK (expr->right.getPtr ());

   bool res = matchAll<N,E,GS> (expr->left, it, sels);
   if (res)  {
      SelSet tmpSels = SelSet::create ();
      tmpSels->append (*sels);
      res = res && matchAll<N,E,GS> (expr->right, it, tmpSels);
   }
   return res;
}

// --- finds one match for SxGQOr
template<class N,class E,template<class,bool> class GS>
bool match (const SxPtr<SxGQOr> &expr,
            const typename SxGraph<N,E,GS>::ConstIterator &it,
            const Selection &sel)
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (expr->left.getPtr ());
   SX_CHECK (expr->right.getPtr ());

   bool res = match<N,E,GS> (expr->left, it, sel);
   if (!res)  {
      res = match<N,E,GS> (expr->right, it, sel);
   }
   return res;

}

// --- finds one match for SxGQOr
template<class N,class E,template<class,bool> class GS>
bool matchAll (const SxPtr<SxGQOr> &expr,
               const typename SxGraph<N,E,GS>::ConstIterator &it,
               const SelSet &sels)
{
   SX_TRACE ();
   SX_CHECK (it.isValid ());
   SX_CHECK (expr->left.getPtr ());
   SX_CHECK (expr->right.getPtr ());

   bool res = matchAll<N,E,GS> (expr->left, it, sels);
   if (!res)  {
      res = matchAll<N,E,GS> (expr->right, it, sels);
   }
   return res;
}

// --- converts to exact expression type
//     and calls corresponding match() function
template<class N,class E,template<class,bool> class GS>
bool match (const SxPtr<SxGQExprBase> &expr,
            const typename SxGraph<N,E,GS>::ConstIterator &it,
            const Selection &sel)
{
   SX_CHECK (expr->isOp () == false); // assuming only AND,OR,==,!=

   SxGQExprBase::ExprType op = expr->exprType;
   if (op == SxGQExprBase::ExprType::AND)  {
      SxPtr<SxGQAnd> p = expr;
      return match<N,E,GS> (p, it, sel);
   }  else if (op == SxGQExprBase::ExprType::OR)  {
      SxPtr<SxGQOr> p = expr;
      return match<N,E,GS> (p, it, sel);
   }  else if (  op == SxGQExprBase::ExprType::Equal
              || op == SxGQExprBase::ExprType::NotEqual
              || op == SxGQExprBase::ExprType::Any)  {
      SxPtr<SxGQExpr> p = expr;
      return match<N,E,GS> (p, it, sel);
   }  else  {
      SX_EXIT; // unsupported type
   }
}
// --- converts to exact expression type
//     and calls corresponding matchAll() function
template<class N,class E,template<class,bool> class GS>
bool matchAll (const SxPtr<SxGQExprBase> &expr,
               const typename SxGraph<N,E,GS>::ConstIterator &it,
               const SelSet &sels)
{
   SX_CHECK (expr->isOp () == false); // assuming only AND,OR,==,!=

   SxGQExprBase::ExprType op = expr->exprType;
   if (op == SxGQExprBase::ExprType::AND)  {
      SxPtr<SxGQAnd> p = expr;
      return matchAll<N,E,GS> (p, it, sels);
   }  else if (op == SxGQExprBase::ExprType::OR)  {
      SxPtr<SxGQOr> p = expr;
      return matchAll<N,E,GS> (p, it, sels);
   }  else if (  op == SxGQExprBase::ExprType::Equal
              || op == SxGQExprBase::ExprType::NotEqual
              || op == SxGQExprBase::ExprType::Any)  {
      SxPtr<SxGQExpr> p = expr;
      return matchAll<N,E,GS> (p, it, sels);
   }  else  {
      SX_EXIT; // unsupported type
   }
}

} // namespace GQ
// -----------------------------------------------------------------------------

// --- returns SxGQExprBase representing an Equal condition
template<class T>
SxPtr<SxGQExprBase>
operator== (const sx::N &n_, const T &val_)
{
   SxPtr<SxGQExpr> e1 = SxPtr<SxGQExpr>::create (n_.getPropName (), SxVariant(val_),
                                                 SxGQExpr::ExprType::Equal);
   return e1;
}

// --- returns SxGQExprBase representing a not-Equal condition
template<class T>
SxPtr<SxGQExprBase>
operator!= (const sx::N &n_, const T &val_)
{
   SxPtr<SxGQExpr> e1 = SxPtr<SxGQExpr>::create (n_.getPropName (), SxVariant(val_),
                                                 SxGQExpr::ExprType::NotEqual);
   return e1;
}

// -----------------------------------------------------------------------------

// --- returns SxGQExprBase representing a binary AND condition
SX_EXPORT_GRAPH SxPtr<SxGQExprBase>
operator&& (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprBase> &e2);

// --- returns SxGQExprBase representing a binary OR condition
SX_EXPORT_GRAPH SxPtr<SxGQExprBase>
operator|| (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprBase> &e2);

// -----------------------------------------------------------------------------
// Ordered Direct

// --- joins two child nodes represented by two SxGQExprBase objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprList> &e2);

// --- joins two child nodes represented by two SxGQExprList objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator>= (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprList> &e2);

// -----------------------------------------------------------------------------
// Unordered Indirect

// --- joins two child nodes represented by two SxGQExpr objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprList> &e2);

// --- joins two child nodes represented by two SxGQExprList objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator< (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprList> &e2);

// -----------------------------------------------------------------------------
// Unordered Direct

// --- joins two child nodes represented by two SxGQExpr objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprList> &e2);

// --- joins two child nodes represented by two SxGQExprList objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator<= (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprList> &e2);

// -----------------------------------------------------------------------------
// Ordered Indirect

// --- joins two child nodes represented by two SxGQExpr objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprList and SxGQExprBase
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprBase> &e2);

// --- joins two child nodes represented by SxGQExprBase and SxGQExprList
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprList> &e2);

// --- joins two child nodes represented by two SxGQExprList objects
SX_EXPORT_GRAPH SxPtr<SxGQExprList>
operator> (const SxPtr<SxGQExprList> &e1, const SxPtr<SxGQExprList> &e2);


// -----------------------------------------------------------------------------

// --- joins SxGQExprBase(parent) with SxGQExprList(children) as SxGQRelations
SX_EXPORT_GRAPH SxPtr<SxGQExprBase>
operator^ (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprList> &elist);

// --- joins SxGQExprBase(parent) with SxGQExprBase(children) as SxGQRelations
SX_EXPORT_GRAPH SxPtr<SxGQExprBase>
operator^ (const SxPtr<SxGQExprBase> &e1, const SxPtr<SxGQExprBase> &e2);

// -----------------------------------------------------------------------------


#endif /*_SX_G_QUERY_H_*/
