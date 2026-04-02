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

#ifndef _SX_GQ_RELATIONS_H_
#define _SX_GQ_RELATIONS_H_

#include <SxGQExprBase.h>
#include <SxGQExprList.h>

/** \brief Graph Query Parent-Child relationship expression

    \b SxGQRelations = SPHInX Graph Query Parent-Child Relationship Class

     SxGQRelations class stores the
     representation for relationship
     between parent and child nodes.

\code

N1 ^ ( N2 < N3 )

where N1 is the parent of N2 and N3

\endcode
 */
class SX_EXPORT_GRAPH SxGQRelations : public SxGQExprBase
{
   public:
      SxGQRelations ();
      SxGQRelations (const SxPtr<SxGQExprBase> &frm_,
                     size_t min_=1, size_t max_=1);
      SxGQRelations (const SxPtr<SxGQExprBase> &frm_,
                     const SxPtr<SxGQExprList> &to_,
                     size_t min_=1, size_t max_=1);
     ~SxGQRelations ();

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

   protected:
      std::ostream &print (std::ostream &s) const override;

      SxPtr<SxGQExprBase> from;
      SxPtr<SxGQExprList> to;
      size_t minDist;
      size_t maxDist;

};

#endif /*_SX_GQ_RELATIONS_H_*/
