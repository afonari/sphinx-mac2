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

#include <SxN.h>

sx::N::N () : propName(""),
              op(SxGQExprBase::ExprType::None),
              captureEnabled(false)
{
   // empty
}

sx::N::N (SxString prop_)
{
   SX_TRACE ();
   SX_CHECK (prop_ != "");
   propName = prop_;
   op = SxGQExprBase::ExprType::None;
   captureEnabled = false;
}

SxString sx::N::getPropName () const
{
   SX_TRACE ();
   return propName;
}

SxString sx::N::getCName () const
{
   SX_TRACE ();
   return captureName;
}

sx::N &sx::N::c (SxString name)
{
   SX_TRACE ();
   captureEnabled = true;
   captureName = name;
   return *this;
}

bool sx::N::isCaptured () const
{
   SX_TRACE ();
   return captureEnabled;
}

sx::N &sx::N::any ()
{
   SX_TRACE ();
   op = SxGQExprBase::ExprType::Any;
   return *this;
}

sx::N::operator SxPtr<SxGQExprBase>() const
{
   SX_TRACE ();
   SxPtr<SxGQExpr> e = SxPtr<SxGQExpr>::create (propName, val,
                                                captureName,
                                                op,
                                                captureEnabled);
   return e;
}


