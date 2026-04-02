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

#ifndef _SX_N_H_
#define _SX_N_H_

#include <SxGQExpr.h>
#include <SxExportGraph.h>

/** \brief Graph Query Node Symbol Class

    \b N = SPHInX Graph Query Node Symbol Class

     N class represents a node symbol
     within a graph query. It's
     constructor takes the property
     name as parameter.

\code
N("prop1") == "val1";

or

N("prop2") != "val2";

or

N("prop3").any();

\endcode
 */
namespace sx {
class SX_EXPORT_GRAPH N
{
   public:
      N ();
      N (SxString prop_);

      SxString getPropName () const;
      SxString getCName () const;

      N &c (SxString name = "");
      bool isCaptured () const;
      N &any ();
      operator SxPtr<SxGQExprBase>() const;

      template<class T>
      N &eq(const T &propVal)
      {
         op  = SxGQExprBase::ExprType::Equal;
         val = SxVariant(propVal);
         return *this;
      }

      template<class T>
      N &neq(const T &propVal)
      {
         op  = SxGQExprBase::ExprType::NotEqual;
         val = SxVariant(propVal);
         return *this;
      }

   protected:
      SxString propName;
      SxVariant val;
      SxString captureName;
      SxGQExprBase::ExprType op;
      bool captureEnabled;
};

}
#endif /*_SX_N_H_*/
