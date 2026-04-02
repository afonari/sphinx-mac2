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

#include <SxGQPattern.h>
#include <SxGQExprBase.h>

size_t sxHash (const SxGQPattern &in)
{
   return SxHashFunction::hash (in.getId ());
}

SxGQPattern::SxGQPattern () : id(-1), relType(RelationType::OrderedDirect)
{
   // empty
}
SxGQPattern::SxGQPattern (ssize_t id_)
   : id(id_), relType(RelationType::OrderedDirect)
{
   // empty
}
SxGQPattern::SxGQPattern (ssize_t id_, const SxPtr<SxGQExprBase> &expr_)
   : id(id_), expr(expr_), relType(RelationType::OrderedDirect)
{
   // empty
}

SxGQPattern::SxGQPattern (const SxGQPattern &in)
{
   SX_TRACE ();
   id      = in.id;
   expr    = in.expr;
   relType = in.relType;
}

SxGQPattern::~SxGQPattern ()
{
   // empty
}

SxGQPattern &SxGQPattern::operator= (const SxGQPattern &in)
{
   SX_TRACE ();
   if (this == &in)  return *this;
   id      = in.id;
   expr    = in.expr;
   relType = in.relType;
   return *this;
}

bool SxGQPattern::operator== (const SxGQPattern &in) const
{
   SX_TRACE ();
   return id == in.id;
}

void SxGQPattern::setRelType (int rType)
{
   SX_TRACE ();
   relType = static_cast<RelationType>(rType);
}

SxGQPattern::RelationType SxGQPattern::getRelType () const
{
   SX_TRACE ();
   return relType;
}

ssize_t SxGQPattern::getId () const
{
   SX_TRACE ();
   return id;
}

void SxGQPattern::setExpr (const SxPtr<SxGQExprBase> &expr_)
{
   SX_TRACE ();
   SX_CHECK (expr_.getPtr ());
   expr = expr_;
}

