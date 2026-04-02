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

#ifndef _SX_GQ_EXPR_H_
#define _SX_GQ_EXPR_H_

#include <SxGQExprBase.h>

/** \brief Graph Query Basic Expression Class

    \b SxGQExpr = SPHInX Graph Query Basic Expression Class

     SxGQExpr class represents a basic
     expression in graph query. In the
     AST this class represents the leaf
     nodes. It matches left and right
     expression for equality or
     non-equality depending on OpType.

\code
#include <SxGQueryNode.h>
#include <SxGQExpr.h>

SxGQExpr expr = (N("propName") == "propVal");

or

SxGQExpr expr = (N("propName") != "propVal");

or

// expr will evaluate to true, if the
// node contains a property with name
// 'propName' with any value.

SxGQExpr expr = (N("propName").any());

\endcode
 */

class SxGQExpr;

// --- functions to evaluate an SxGQExpr
namespace SxGQInternal {

   template<class N,class E,template<class,bool> class GS>
   bool match (const SxPtr<SxGQExpr> &expr,
               const typename SxGraph<N,E,GS>::ConstIterator &it,
               const SxGQExprBase::Selection &sel);

   template<class N,class E,template<class,bool> class GS>
   bool matchAll (const SxPtr<SxGQExpr> &expr,
                  const typename SxGraph<N,E,GS>::ConstIterator &it,
                  const SxGQExprBase::SelSet &sels);
}

class SX_EXPORT_GRAPH SxGQExpr : public SxGQExprBase
{
   public:
      SxGQExpr ();
      SxGQExpr (const SxString &propName_,
                const SxVariant &right_,
                const SxGQExprBase::ExprType &type_);
      SxGQExpr (const SxString &propName_,
                const SxVariant &right_,
                const SxString &capName_,
                const SxGQExprBase::ExprType &type_,
                bool isCap_);


      SxPtr<SxList<SxPtr<SxGQExprBase> > > firsts () const override;
      SxPtr<SxList<SxPtr<SxGQExprBase> > > lasts () const override;
      SxPtr<SxGQExprBase> first() const override;
      SxPtr<SxGQExprBase> last() const override;

      void setLast (const SxPtr<SxGQExprBase> &p) override;
      void setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst) override;

      bool isOp () const override;
      SxGQExprBase::ExprType getRightOp () const override;

      bool isCaptured () const;
      size_t getHash () const override;

      void makeGraph (SxGraph<SxGQPattern> *gPtr) const override;
      SxGQPattern getGraphNode () const override;

      template<class N,class E,template<class,bool> class GS>
      friend bool SxGQInternal::match (const SxPtr<SxGQExpr> &expr,
                                       const typename SxGraph<N,E,GS>::ConstIterator &it,
                                       const Selection &sel);

      template<class N,class E,template<class,bool> class GS>
      friend bool SxGQInternal::matchAll (const SxPtr<SxGQExpr> &expr,
                                          const typename SxGraph<N,E,GS>::ConstIterator &it,
                                          const SelSet &sels);

   protected:

      std::ostream &print (std::ostream &s) const override;

      SxString propName;
      SxString captureName;
      bool captured;
      SxVariant right;
      ssize_t id;

};


#endif /*_SX_GQ_EXPR_H_*/
