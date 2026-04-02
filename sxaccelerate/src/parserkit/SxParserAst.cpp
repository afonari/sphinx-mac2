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

#include <SxParserAst.h>
#include <SxSed.h>

SxParserAst::ElemType
SxParserAst::getTypeId (const SxString &str)
{
   if (str == "int")  {
      return ElemType::Int;
   }  else if (str == "real")  {
      return ElemType::Double;
   }  else if (str == "bool")  {
      return ElemType::Bool;
   }  else if (str == "string")  {
      return ElemType::String;
   }  else if (str == "list")  {
      return ElemType::List;
   }  else if (str == "group")  {
      return ElemType::Group;
   }  else if (str == "vector")  {
      return ElemType::Vector;
   }  else if (str == "matrix")  {
      return ElemType::Matrix;
   }  else if (str == "enum")  {
      return ElemType::Enum;
   }

   return ElemType::Undefined;
}

SxString SxParserAst::getTypeStr (const SxParserAst::ElemType &type_)
{
   switch(type_)  {
      case ElemType::Int:
         return "int";
      case ElemType::Double:
         return "real";
      case ElemType::Bool:
         return "bool";
      case ElemType::String:
         return "string";
      case ElemType::List:
         return "list";
      case ElemType::Group:
         return "group";
      case ElemType::Vector:
         return "vector";
      case ElemType::Matrix:
         return "matrix";
      case ElemType::Enum:
         return "enum";
      case ElemType::Undefined:
         return "undefined";
      default:
         return "unknown";
   }
}

SxString SxParserAst::escapeStr (const SxString &str)
{
   SxSed sedQuote("(\")", "\\\\$1", "g");
   SxSed sedSlash("(\\\\)", "\\\\$1", "g");
   SxSed sedSlashT("(\\t)", "\\\\t", "g");
   SxSed sedSlashR("(\\r)", "\\\\r", "g");
   SxSed sedSlashA("(\\a)", "\\\\a", "g");
   SxSed sedSlashF("(\\f)", "\\\\f", "g");

   SxString res = str;
   res = sedSlash.subst (res);
   res = res.substitute ("\n","\\n");
   res = sedSlashT.subst (res);
   res = sedSlashR.subst (res);
   res = sedSlashA.subst (res);
   res = sedSlashF.subst (res);
   res = sedQuote.subst (res);

   return res;
}

SxString SxParserAst::unescapeStr (const SxString &str)
{
   SxSed sedQuote("(\\\\\")", "\"", "g");
   SxSed sedSlash("(\\\\\\\\)", "\\", "g");
   SxSed sedSlashN("(?<!\\\\)((\\\\\\\\)*)(\\\\n)", "$1\\n", "g");
   SxSed sedSlashT("(?<!\\\\)((\\\\\\\\)*)(\\\\t)", "$1\\t", "g");
   SxSed sedSlashR("(?<!\\\\)((\\\\\\\\)*)(\\\\r)", "$1\\r", "g");
   SxSed sedSlashA("(?<!\\\\)((\\\\\\\\)*)(\\\\a)", "$1\\a", "g");
   SxSed sedSlashF("(?<!\\\\)((\\\\\\\\)*)(\\\\f)", "$1\\f", "g");


   SxString res = str;
   res = sedQuote.subst (res);
   res = sedSlashF.subst (res);
   res = sedSlashA.subst (res);
   res = sedSlashR.subst (res);
   res = sedSlashT.subst (res);
   res = sedSlashN.subst (res);
   res = sedSlash.subst (res);

   return res;
}

// ---------------------------------------------------------------------------
SxParserAst::SxParserAst () : idx(-1)
{
   SX_TRACE ();
   ast = SxPtr<SxGraph<SxGProps> >::create ();
}

SxParserAst::~SxParserAst ()
{
   SX_TRACE ();
}

SxGProps *SxParserAst::addNode (int type_)
{
   SX_TRACE ();
   SxGProps n(++idx);
   n.setProperty (".type", type_);
   n.setProperty (".key", "");
   n.setProperty (".val", SxVariant ());
   return const_cast<SxGProps *>(&(*(ast->createNode (n))));
}

void SxParserAst::addEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast->createEdge (SxGProps(a), SxGProps(b));
}

void SxParserAst::addEdge (const SxGProps *a, const SxGProps *b)
{
   SX_TRACE ();
   ast->createEdge (*a, *b);
}

void SxParserAst::prependEdge (ssize_t a, ssize_t b)
{
   SX_TRACE ();
   ast->createEdge (SxGProps(a), SxGProps(b));
}

SxGProps *SxParserAst::getNode (ssize_t i)
{
   SX_TRACE ();
   return &(*(ast->begin (sx::Forward, SxGProps(i))));
}

const SxPtr<SxGraph<SxGProps> > &SxParserAst::getAst () const
{
   return ast;
}

void SxParserAst::printAst () const
{
   SX_TRACE ();
   ast->print (ast->begin(), cout);
}

