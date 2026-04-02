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

#ifndef _SX_GQ_AND_H_
#define _SX_GQ_AND_H_

#include <SxGQExprBase.h>

/** \brief Graph Query Binary And operator

    \b SxGQAnd = SPHInX Graph Query Binary And Operator Class

     SxGQAnd class represents a binary AND
     operator between two nodes in a query.
     A binary AND evaluatese to true when
     both left and right operands evaluate
     to true.

\code
N1 && N2 && N3 ...

where N is object of class sx::N
\endcode
 */

class SxGQAnd;

// --- functions to evaluate SxGQAnd
namespace SxGQInternal
{
   template<class N,class E,template<class,bool> class GS>
   bool match (const SxPtr<SxGQAnd> &expr,
               const typename SxGraph<N,E,GS>::ConstIterator &it,
               const SxGQExprBase::Selection &sel);

   template<class N,class E,template<class,bool> class GS>
   bool matchAll (const SxPtr<SxGQAnd> &expr,
                  const typename SxGraph<N,E,GS>::ConstIterator &it,
                  const SxGQExprBase::SelSet &sels);
}

class SX_EXPORT_GRAPH SxGQAnd : public SxGQExprBase
{
   public:
      SxGQAnd ();
      SxGQAnd (const SxPtr<SxGQExprBase> &l, const SxPtr<SxGQExprBase> &r);
     ~SxGQAnd ();


      SxPtr<SxList<SxPtr<SxGQExprBase> > > firsts () const override;
      SxPtr<SxList<SxPtr<SxGQExprBase> > > lasts () const override;
      SxPtr<SxGQExprBase> first () const override;
      SxPtr<SxGQExprBase> last () const override;
      void setLast (const SxPtr<SxGQExprBase> &p) override;
      void setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &lst) override;

      bool isOp () const override;
      SxGQExprBase::ExprType getRightOp () const override;
      size_t getHash () const override;

      void makeGraph (SxGraph<SxGQPattern> *gPtr) const override;
      SxGQPattern getGraphNode () const override;

      template<class N,class E,template<class,bool> class GS>
      friend bool SxGQInternal::match (const SxPtr<SxGQAnd> &expr,
                                       const typename SxGraph<N,E,GS>::ConstIterator &it,
                                       const Selection &sel);

      template<class N,class E,template<class,bool> class GS>
      friend bool SxGQInternal::matchAll (const SxPtr<SxGQAnd> &expr,
                                          const typename SxGraph<N,E,GS>::ConstIterator &it,
                                          const SelSet &sels);


   protected:
      std::ostream &print (std::ostream &s) const override;

      SxPtr<SxGQExprBase> left;
      SxPtr<SxGQExprBase> right;
};


#endif /* _SX_GQ_AND_H_ */
