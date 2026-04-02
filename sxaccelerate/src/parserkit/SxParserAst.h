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

#ifndef _SX_PARSER_AST_H_
#define _SX_PARSER_AST_H_

#include <SxParserKit.h>
#include <SxString.h>
#include <SxGraph.h>
#include <SxGProps.h>
/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author John Doe, johndoe@example.org */
class SX_EXPORT_PARSER_KIT SxParserAst
{
   public:

      enum ElemType { Undefined = SxVariantType::Undefined,
                      Int = SxVariantType::Int,
                      Double = SxVariantType::Double,
                      Bool = SxVariantType::Bool,
                      String = SxVariantType::String,
                      List = SxVariantType::List, Group,
                      Vector, Matrix, Enum };

      SxParserAst ();
     ~SxParserAst ();

      SxGProps *addNode (int type_);

      void addEdge (ssize_t, ssize_t);
      void addEdge (const SxGProps *, const SxGProps *);
      void prependEdge (ssize_t, ssize_t);

      SxGProps *getNode (ssize_t);

      const SxPtr<SxGraph<SxGProps> > &getAst () const;

      void printAst () const;

      static ElemType getTypeId (const SxString &str);
      static SxString getTypeStr (const ElemType &type_);

      static SxString escapeStr (const SxString &str);
      static SxString unescapeStr (const SxString &str);

   protected:

      ssize_t idx;
      SxPtr<SxGraph<SxGProps> > ast;

};

#endif /* _SX_PARSER_AST_H_ */
