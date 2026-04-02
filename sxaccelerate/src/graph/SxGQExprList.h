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

#ifndef _SX_GQ_EXPR_LIST_H_
#define _SX_GQ_EXPR_LIST_H_

#include <SxGQExprBase.h>

/** \brief Graph Query List of expressions

    \b SxGQExprList = SPHInX Graph Query Expression List Class

    SxGQExprList class stores list of expressions
    which can be of type SxGQExprList, SxGQExpr and
    SxGQRelations.

\code
(N1 op N2 op N3 ...)

where op can be either of the four
types of siblings: >=,<,<=,>
\endcode
 */
class SX_EXPORT_GRAPH SxGQExprList : public SxGQExprBase
{

   public:
      enum SiblingType {OrderedDirect, OrderedIndirect,
                        UnorderedDirect, UnorderedIndirect};
      SxGQExprList ();

      SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                    SiblingType sType = OrderedDirect);

      SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                    const SxPtr<SxGQExprBase> &e2,
                    SiblingType sType = OrderedDirect);

      SxGQExprList (const SxPtr<SxGQExprBase> &e1,
                    const SxPtr<SxGQExprList> &e2,
                    SiblingType sType = OrderedDirect);

      SxGQExprList (const SxPtr<SxGQExprList> &e1,
                    const SxPtr<SxGQExprBase> &e2,
                    SiblingType sType = OrderedDirect);

     ~SxGQExprList ();

      void append (const SxGQExprList &lst);
      void append (const SxPtr<SxGQExprBase> &e);
      void append (const SxPtr<SxGQExprList> &lst);

      SiblingType getSibType () const;
      const SxList<SxPtr<SxGQExprBase> > &getExprList () const;


      SxPtr<SxList<SxPtr<SxGQExprBase> > > firsts () const override;
      SxPtr<SxList<SxPtr<SxGQExprBase> > > lasts () const override;
      void setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst) override;
      SxPtr<SxGQExprBase> first () const override;
      SxPtr<SxGQExprBase> last () const override;
      void setLast (const SxPtr<SxGQExprBase> &p) override;

      bool isOp () const override;
      SxGQExprBase::ExprType getRightOp () const override;

      size_t getHash () const override;

      void makeGraph (SxGraph<SxGQPattern> *gPtr) const override;
      SxGQPattern getGraphNode () const override;

   protected:
      std::ostream &print (std::ostream &s) const override;

      SxList<SxPtr<SxGQExprBase>> exprList;
      const SiblingType sibType;
};

#endif /*_SX_GQ_EXPR_LIST_H_*/
