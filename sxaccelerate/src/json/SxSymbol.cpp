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

#include <SxSymbol.h>
#include <SxFileIO.h>
#include <SxUUIDv4.h>

namespace SxParserKit {

SxSymbol::SxSymbol ()
   : name (""),
     type (Type::Undefined),
     dataG (),
     nodeId (-1)
{
   dataG = SxPtr<SxGraph<SxGProps> >::create ();
   nodeId = dataG->createNode (SxGProps (getIdx ())).getIdx ();
   auto it = dataG->begin (nodeId);
   it->setProperty (".type", (int)Type::Undefined);
   it->setProperty (".key", "");
   it->setProperty (".val", SxVariant ());
}

SxSymbol::SxSymbol (const SxPtr<SxGraph<SxGProps> > &gPtr_,
                    ssize_t nodeId_)
{
   SX_CHECK (gPtr_.getPtr ());

   dataG   = gPtr_;
   nodeId  = nodeId_;

   auto it = dataG->begin (nodeId);

   SX_CHECK (it.isValid ());
   if (it.isValid ()) {
      name  = it->getProperty (".key").toString ();
      type  = (Type)it->getProperty (".type").getInt ();
   } else {
      name    = "";
      type    = Type::Undefined;
      nodeId_ = -1;
   }
}

SxSymbol::~SxSymbol ()
{
   // empty
}

SxString SxSymbol::toString () const
{
   SX_TRACE ();

   if (type != Type::String)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Element of type '")
                + SxParserAst::getTypeStr (type)
                + "' cannot be converted to 'string'");

   return dataG->begin (nodeId)->getProperty (".val").toString ();
}

int64_t SxSymbol::toInt () const
{
   SX_TRACE ();

   bool isAllowed = (type == Type::Int || type == Type::Double);

   if (!isAllowed)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Element '")
                + name + "' of type '"
                + SxParserAst::getTypeStr (type)
                + "' cannot be converted to 'int'");

   return dataG->begin (nodeId)->getProperty (".val").toInt ();
}

double SxSymbol::toDouble () const
{
   SX_TRACE ();

   bool isAllowed = (type == Type::Int || type == Type::Double);

   if (!isAllowed)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Element '")
                + name + "' of type '"
                + SxParserAst::getTypeStr (type)
                + "' cannot be converted to 'double'");

   return dataG->begin (nodeId)->getProperty (".val").toDouble ();
}

bool SxSymbol::toBool () const
{
   SX_TRACE ();

   if (type != Type::Bool)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Element '")
                + name + "' of type '"
                + SxParserAst::getTypeStr (type)
                + "' cannot be converted to 'bool'");

   return dataG->begin (nodeId)->getProperty (".val").toBool ();
}

SxList<int64_t> SxSymbol::toIntList () const
{
   SX_TRACE ();

   if (type != Type::List)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'list', got '")
                + SxParserAst::getTypeStr (type) + "'");

   SxList<int64_t> res;
   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      int nType = (int)chIt->getProperty (".type").getInt ();
      if (nType == Type::Double || nType == Type::Int)  {
         res.append (chIt->getProperty (".val").toInt ());
      }
   }
   return res;
}

SxList<SxString> SxSymbol::toStringList () const
{
   SX_TRACE ();

   if (type != Type::List)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'list', got '")
                + SxParserAst::getTypeStr (type) + "'");

   SxList<SxString> res;
   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      int nType = (int)chIt->getProperty (".type").getInt ();
      if (nType == Type::String)  {
         res.append (chIt->getProperty (".val").toString ());
      }
   }
   return res;
}

SxList<double> SxSymbol::toDoubleList () const
{
   SX_TRACE ();

   if (type != Type::List)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'list', got '")
                + SxParserAst::getTypeStr (type) + "'");

   SxList<double> res;
   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      int nType = (int)chIt->getProperty (".type").getInt ();
      if (nType == Type::Double || nType == Type::Int)  {
         res.append (chIt->getProperty (".val").toDouble ());
      }
   }
   return res;
}

SxList<SxSymbol> SxSymbol::toList () const
{
   SX_TRACE ();

   if (type != Type::List && type != Type::Group)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'list' or 'group', got '")
                + SxParserAst::getTypeStr (type) + "'");

   SxList<SxSymbol> res;
   // --- collect all child elements
   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      res.append (SxSymbol (dataG, chIt.getIdx ()));
   }

   return res;
}

bool SxSymbol::hasElem (const SxString &key) const
{
   SX_TRACE ();

   if (type != Type::Group)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'group', got '")
                + SxParserAst::getTypeStr (type) + "'");

   auto it = dataG->begin (nodeId);
   // --- search within all child elements
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      if (chIt->getProperty (".key").toString () == key)  {
         return true;
      }
   }

   return false;
}

SxSymbol SxSymbol::getElem (const SxString &key, bool isRequired) const
{
   SX_TRACE ();

   if (type != Type::Group)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'group', got '")
                + SxParserAst::getTypeStr (type) + "'");

   if (isRequired && !hasElem (key))
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Element with key '") + key
                + "' not found in group '" + name + "'");

   auto it = dataG->begin (nodeId);
   // --- search within all child elements
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      if (chIt->getProperty (".key").toString () == key)  {
         SxSymbol s(dataG, chIt.getIdx ());
         return s;
      }
   }
   return SxSymbol ();
}

SxSymbol SxSymbol::getElem (ssize_t idx) const
{
   if (type != Type::List)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'list', got '")
               + SxParserAst::getTypeStr (type) + "'");

   auto it = dataG->begin (nodeId);
   if (idx < 0 || idx >= it.getSizeOut ()) {
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("List index out of range ")
               + "0 <= " + idx + " < " + it.getSizeOut ());
   }

   return SxSymbol (dataG, it.out (idx).getIdx ());
}

SxList<SxSymbol> SxSymbol::getTypedList (Type type_, bool recursive) const
{
   SX_TRACE ();

   if (type != Type::List && type != Type::Group)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'list' or 'group', got '")
                + SxParserAst::getTypeStr (type) + "'");

   SxList<SxSymbol> res;

   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      SxSymbol tr (dataG, chIt.getIdx ());
      if (tr.getType () == type_)  {
         res.append (SxSymbol (dataG, chIt.getIdx ()));
         if (recursive)  {
            Type t = tr.getType ();
            if (t == Type::Group || t == Type::List)  {
               res << tr.getTypedList (type_, recursive);
            }
         }
      }
   }
   return res;
}

SxList<SxSymbol>
SxSymbol::getTypedList (const SxString &key_, Type type_, bool recursive) const
{
   SX_TRACE ();

   if (type != Type::Group)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'group', got '")
                + SxParserAst::getTypeStr (type) + "'");

   SxList<SxSymbol> res;

   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i)  {
      auto chIt = it.out (i);
      SxSymbol tr (dataG, chIt.getIdx ());
      if (tr.getType () == type_ && tr.getKey () == key_)  {
         res.append (SxSymbol (dataG, chIt.getIdx ()));
         if (recursive)  {
            Type t = tr.getType ();
            if (t == Type::Group || t == Type::List)  {
               res << tr.getTypedList (key_, type_, recursive);
            }
         }
      }
   }

   return res;
}

SxList<SxSymbol> SxSymbol::getGroupList (bool recursive) const
{
   SX_TRACE ();
   return getTypedList (Type::Group, recursive);
}

SxList<SxSymbol> SxSymbol::getGroupList (const SxString &key_,
                                         bool recursive) const
{
   SX_TRACE ();
   return getTypedList (key_, Type::Group, recursive);
}

SxList<SxSymbol> SxSymbol::getArrayList (bool recursive) const
{
   SX_TRACE ();
   return getTypedList (Type::List, recursive);
}

SxList<SxSymbol> SxSymbol::getArrayList (const SxString &key_,
                                         bool recursive) const
{
   SX_TRACE ();
   return getTypedList (key_, Type::List, recursive);
}

SxString SxSymbol::getKey () const
{
   SX_TRACE ();
   return name;
}

SxVariant SxSymbol::getValue () const
{
   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));

   if (type == Type::Group || type == Type::List)
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not of primitive type"));

   return dataG->begin (nodeId)->getProperty (".val");
}

SxSymbol::Type SxSymbol::getType () const
{
   SX_TRACE ();
   return type;
}

SxString SxSymbol::getTag () const
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));
   return dataG->begin (nodeId)->getProperty (".val").getTag ();
}

SxString SxSymbol::getDocTxt () const
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));
   if (dataG->begin (nodeId)->hasProperty (".docTxt"))
      return dataG->begin (nodeId)->getProperty (".docTxt").getString ();
   return "";
}

SxString SxSymbol::getDocTxtTag () const
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));

   auto it = dataG->begin (nodeId);

   if (it->hasProperty (".docTxt"))
      return it->getProperty (".docTxt").getTag ();
   return "";
}

bool SxSymbol::isValid () const
{
   SX_TRACE ();
   if (type == Type::Undefined)  return false;
   // in case node has been deleted
   return (dataG->begin (nodeId).isValid ());
}

void SxSymbol::print (ostream &os, const SxString &newLine,
                      int lvl, int tabs,
                      bool isFileStream) const
{
   SX_TRACE ();
   SX_CHECK (isValid ());
   SX_UNUSED (newLine);

   if (type == Type::Group)  {
      printObject (os, newLine, lvl, tabs, isFileStream);
   }  else if (type == Type::List)  {
      printArray (os, newLine, lvl, tabs, isFileStream);
   }  else if (type == Type::Int)  {
      printInt (os, newLine, lvl, tabs, isFileStream);
   }  else if (type == Type::Double)  {
      printDouble (os, newLine, lvl, tabs, isFileStream);
   }  else if (type == Type::String)  {
      printString (os, newLine, lvl, tabs, isFileStream);
   }  else if (type == Type::Bool)  {
      printBool (os, newLine, lvl, tabs, isFileStream);
   }  else  {
      SX_EXIT;
   }
}

void SxSymbol::printObject (ostream &os, const SxString &newLine,
                            int lvl, int tabs,
                            bool isFileStream) const
{
   SX_TRACE ();
   SX_CHECK (type == Type::Group, type);
   SxList<SxSymbol> elems = toList ();
   ssize_t c = 0;

   if (getDocTxt () != "") os << getDocTxt () << newLine;

   for (int i = 0; i < lvl*tabs; ++i)  os << " ";

   if (name != "" && !isRoot ())  {
      if (isFileStream)  os << "\"" << SxParserAst::escapeStr (name) << "\":";
      else               os << "\"" << name                          << "\":";
   }
   os << "{";
   for (auto it = elems.begin (); it != elems.end (); ++it)
   {
      if (c != 0)  os << ",";
      os << newLine;
      it->print (os, newLine, lvl+1, tabs, isFileStream);
      c++;
   }
   os << newLine;
   for (int i = 0; i < lvl*tabs; ++i)  os << " ";
   os << "}";
}

void SxSymbol::printArray (ostream &os, const SxString &newLine,
                           int lvl, int tabs,
                           bool isFileStream) const
{
   SX_TRACE ();
   SX_CHECK (type == Type::List, type);
   SxList<SxSymbol> elems = toList ();
   ssize_t c = 0;

   if (getDocTxt () != "") os << getDocTxt () << newLine;

   for (int i = 0; i < lvl*tabs; ++i)  os << " ";

   if (name != "" && !isRoot ())  {
      if (isFileStream)  os << "\"" << SxParserAst::escapeStr (name) << "\":";
      else               os << "\"" << name                          << "\":";
   }

   os << "[";
   for (auto it = elems.begin (); it != elems.end (); ++it)
   {
      if (c != 0)  os << ",";
      os << newLine;
      it->print (os, newLine, lvl+1, tabs, isFileStream);
      c++;
   }
   os << newLine;
   for (int i = 0; i < lvl*tabs; ++i)  os << " ";
   os << "]";
}

void SxSymbol::printInt (ostream &os, const SxString &newLine,
                         int lvl, int tabs,
                         bool isFileStream) const
{
   SX_TRACE ();
   SX_CHECK (type == Type::Int, type);

   if (getDocTxt () != "") os << getDocTxt () << newLine;

   for (int i = 0; i < lvl*tabs; ++i)  os << " ";

   if (name != "")  {
      if (isFileStream)  os << "\"" << SxParserAst::escapeStr (name) << "\":";
      else               os << "\"" << name                          << "\":";
   }

   os << toInt ();
}

void SxSymbol::printDouble (ostream &os, const SxString &newLine,
                            int lvl, int tabs,
                            bool isFileStream) const
{
   SX_TRACE ();
   SX_CHECK (type == Type::Double, type);

   if (getDocTxt () != "") os << getDocTxt () << newLine;

   for (int i = 0; i < lvl*tabs; ++i)  os << " ";

   if (name != "")  {
      if (isFileStream)  os << "\"" << SxParserAst::escapeStr (name) << "\":";
      else               os << "\"" << name                          << "\":";
   }

   os << toDouble ();
}

void SxSymbol::printString (ostream &os, const SxString &newLine,
                            int lvl, int tabs,
                            bool isFileStream) const
{
   SX_TRACE ();
   SX_CHECK (type == Type::String, type);

   if (getDocTxt () != "") os << getDocTxt () << newLine;

   for (int i = 0; i < lvl*tabs; ++i)  os << " ";

   if (isFileStream)  {
      if (name != "")  os << "\"" << SxParserAst::escapeStr (name) << "\":";
      os << "\"" << SxParserAst::escapeStr (toString ()) << "\"";
   }  else  {
      if (name != "")  os << "\"" << name << "\":";
      os << "\"" << toString () << "\"";
   }
}

void SxSymbol::printBool (ostream &os, const SxString &newLine,
                          int lvl, int tabs,
                          bool isFileStream) const
{
   SX_TRACE ();
   SX_CHECK (type == Type::Bool, type);

   if (getDocTxt () != "") os << getDocTxt () << newLine;

   for (int i = 0; i < lvl*tabs; ++i)  os << " ";

   if (name != "")  {
      if (isFileStream)  os << "\"" << SxParserAst::escapeStr (name) << "\":";
      else               os << "\"" << name                          << "\":";
   }

   if (toBool ())  os << "true";
   else            os << "false";
}

void SxSymbol::setElem (ssize_t idx, const SxVariant &val)
{
   SX_TRACE ();

   if (type != Type::List)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'list', got '")
                + SxParserAst::getTypeStr (type) + "'");

   if (val.getType () == Type::List)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("SxVariant of type 'list' cannot be set"));

   auto it = dataG->begin (nodeId);
   if (idx < 0 || idx >= it.getSizeOut ()) {
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("List index out of range ")
               + "0 <= " + idx + " < " + it.getSizeOut ());
   }
   SxSymbol sym (dataG, it.out (idx).getIdx ());
   sym.setValue (val);
}

void SxSymbol::setElem (const SxString &key, const SxVariant &val)
{
   SX_TRACE ();

   if (type != Type::Group)
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Expected type 'group', got '")
                + SxParserAst::getTypeStr (type) + "'");

   if (!hasElem (key))
      SX_THROW (  "SxSymbol", "InvalidOperation",
                  SxString ("Key '") + key + "' does not exists");

   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i) {
      auto elemIt = it.out (i);
      if (elemIt->getProperty (".key").getString () == key) {
         SxSymbol sym(dataG, elemIt.getIdx ());
         sym.setValue (val);
         break;
      }
   }
}

void SxSymbol::setType (Type type_)
{
   SX_TRACE ();

   if (type_ == type)  return;

   if (type_ != Type::Group && type_ != Type::List) {
      // root has to be group or a list
      if (isRoot () == true) {
         SX_THROW ("SxSymbol", "InvalidOperation",
                   "root element must be array or group");
      }
   }

   auto it = dataG->begin (nodeId);
   if (type == Type::List || type == Type::Group) {
      // remove the child elements
      SxPtr<SxUniqueList<ssize_t> > selection =
         SxPtr<SxUniqueList<ssize_t> >::create ();
      for (ssize_t i = 0; i < it.getSizeOut (); ++i) {
         selection->append (it.out (i).getIdx ());
      }
      auto selIt = selection->begin ();
      for (; selIt.isValid (); ++selIt) {
         removeRecursive (*selIt);
      }
   }

   type = type_;
   it->setProperty (".type", (int)type);
   it->setProperty (".val", SxVariant());
}

void SxSymbol::setKey (const SxString &key)
{
   SX_TRACE ();

   if (key == "")
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Key cannot be an empty string"));

   // check if parent exists and is a group
   auto it = dataG->begin (nodeId);
   if (it.getSizeIn () > 0) {
      auto parentIt = it.in (0);
      if (parentIt->getProperty (".type").getInt () == (int)Type::Group) {
         it->setProperty (".key", key);
      } else {
         SX_THROW ("SxSymbol", "IOError",
                   "A key is only allowed when parent is a group");
      }
   } else {
      SX_THROW ("SxSymbol", "IOError",
                "Root element cannot have a key");
   }
}

void SxSymbol::setValue (const SxVariant &val)
{
   SX_TRACE ();

   if (!val.isInitialized ())
      SX_THROW ("SxSymbol", "InvalidInput",
                "uninitialized value cannot be set");

   setType ((Type)val.getType ());
   dataG->begin (nodeId)->setProperty (".val", val);
}

void SxSymbol::setDocTxt (const SxString &str, const SxString &tag)
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));

   if (str != "" && !isValidDocTxt (str))
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Given doc text does not match required pattern"));

   SxVariant v(str);
   v.setTag (tag);
   dataG->begin (nodeId)->setProperty (".docTxt", v);
}

SxSymbol SxSymbol::append (Type type_)
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Cannot append to uninitialized symbol"));

   if (type_ == Type::Undefined)
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Cannot append an element of type Undefined"));

   if (type != Type::List)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'list', got '")
               + SxParserAst::getTypeStr (type) + "'");

   SxGProps n (getIdx ());
   n.setProperty (".type", (int)type_);
   n.setProperty (".key", "");
   n.setProperty (".val", SxVariant());
   auto cIt = dataG->createNode (n);
   SX_CHECK (cIt.isValid ());
   dataG->createEdge (*(dataG->begin(nodeId)), *cIt);
   return SxSymbol (dataG, cIt.getIdx ());
}

SxSymbol SxSymbol::append (const SxString &key, Type type_)
{
   SX_TRACE ();

   if (key == "")
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Key cannot be an empty string"));

   if (type_ == Type::Undefined)
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Cannot append an element of type Undefined"));

   if (type != Type::Group)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'group', got '")
               + SxParserAst::getTypeStr (type) + "'");

   if (hasElem (key))
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Key '") + key + "' already exists");

   SxGProps n (getIdx ());
   n.setProperty (".type", (int)type_);
   n.setProperty (".key", key);
   n.setProperty (".val", SxVariant());
   auto cIt = dataG->createNode (n);
   SX_CHECK (cIt.isValid ());
   dataG->createEdge (*(dataG->begin(nodeId)), *cIt);
   return SxSymbol (dataG, cIt.getIdx ());
}

SxSymbol SxSymbol::append (const SxVariant &val)
{
   SX_TRACE ();

   if (type != Type::List)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'list', got '")
               + SxParserAst::getTypeStr (type) + "'");

   if (!val.isInitialized ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                "Cannot append uninitialized value");

   // cannot append a SxVariant(list) directly
   // use append(Type::List)
   if (val.getType () == Type::List)
      SX_THROW ( "SxSymbol", "InvalidInput",
                 SxString ("SxVariant of type 'list' cannot be appended"));

   SxGProps n (getIdx ());
   n.setProperty (".type", (int)val.getType ());
   n.setProperty (".key", "");
   n.setProperty (".val", val);
   auto cIt = dataG->createNode (n);
   SX_CHECK (cIt.isValid ());
   dataG->createEdge (*(dataG->begin(nodeId)), *cIt);
   return SxSymbol (dataG, cIt.getIdx ());
}

SxSymbol SxSymbol::append (const SxString &key, const SxVariant &val)
{
   SX_TRACE ();

   if (key == "")
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Key cannot be an empty string"));

   if (!val.isInitialized ())
      SX_THROW ("SxSymbol", "InvalidInput",
                "Cannot append uninitialized value");

   if (type != Type::Group)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'group', got '")
               + SxParserAst::getTypeStr (type) + "'");

   if (hasElem (key))
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Key '") + key + "' already exists");

   SxGProps n (getIdx ());
   n.setProperty (".type", (int)val.getType ());
   n.setProperty (".key", key);
   n.setProperty (".val", val);
   auto cIt = dataG->createNode (n);
   SX_CHECK (cIt.isValid ());
   dataG->createEdge (*(dataG->begin(nodeId)), *cIt);
   return SxSymbol (dataG, cIt.getIdx ());
}

void SxSymbol::removeRecursive (ssize_t nodeId_)
{
   SX_TRACE ();
   SxPtr<SxUniqueList<ssize_t> > selection =
      SxPtr<SxUniqueList<ssize_t> >::create ();
   auto it = dataG->begin (nodeId_, sx::Forward);
   while (it.isValid ()) {
      selection->append (it.getIdx ());
      ++it;
   }

   auto selIt = selection->begin ();
   for (; selIt.isValid (); ++selIt) {
      dataG->removeNode (*selIt);
   }
}

void SxSymbol::remove (const SxString &key)
{
   SX_TRACE ();

   if (type != Type::Group)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'group', got '")
               + SxParserAst::getTypeStr (type) + "'");

   if (key == "")
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Key cannot be an empty string"));

   if (!hasElem (key))
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Key '") + key + "' does not exist");

   auto it = dataG->begin (nodeId);
   for (ssize_t i = 0; i < it.getSizeOut (); ++i) {
      auto elemIt = it.out (i);
      if (elemIt->getProperty (".key").getString () == key) {
         removeRecursive (elemIt.getIdx ());
         break;
      }
   }
}

void SxSymbol::remove (ssize_t idx)
{
   SX_TRACE ();

   if (type != Type::List)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'list', got '")
               + SxParserAst::getTypeStr (type) + "'");

   auto it = dataG->begin (nodeId);
   if (idx < 0 || idx >= it.getSizeOut ()) {
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("List index out of range ")
               + "0 <= " + idx + " < " + it.getSizeOut ());
   }
   removeRecursive (it.out (idx).getIdx ());
}

SxSymbol SxSymbol::find (const SxString &key) const
{
   SX_TRACE ();

   if (key == "")
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Key cannot be an empty string"));

   if (type != Type::List && type != Type::Group)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'list' or 'group', got '")
               + SxParserAst::getTypeStr (type) + "'");

   SxGQuery q = (sx::N(".key") == key);
   SxGQuery::Selection sel = q.match (dataG, dataG->begin (nodeId));
   if (sel->getSize () > 0) {
      return SxSymbol (dataG, (*sel)(0));
   }
   return SxSymbol ();
}

SxList<SxSymbol> SxSymbol::findAll (const SxString &key) const
{
   SX_TRACE ();

   if (key == "")
      SX_THROW ("SxSymbol", "InvalidInput",
                SxString ("Key cannot be an empty string"));

   if (type != Type::List && type != Type::Group)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'list' or 'group', got '")
               + SxParserAst::getTypeStr (type) + "'");

   SxList<SxSymbol> res;
   SxGQuery q = (sx::N(".key") == key);
   SxGQuery::SelSet sels = q.matchAll (dataG, dataG->begin (nodeId));
   for (auto selsIt = sels->begin (); selsIt.isValid (); ++selsIt) {
      for (auto it2 = dataG->begin (*selsIt); it2.isValid (); ++it2) {
         res.append (SxSymbol (dataG, it2.getIdx ()));
      }
   }
   return res;
}

SxSymbol SxSymbol::getParent () const
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));

   if (isRoot ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Root element cannot have a parent"));

   auto it = dataG->begin (nodeId);
   if (it.getSizeIn () > 0) {
      return SxSymbol (dataG, it.in (0).getIdx ());
   }
   return SxSymbol ();
}

SxList<SxString> SxSymbol::getPath () const
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));

   SxList<SxString> res;
   auto it = dataG->begin (nodeId).reverse ();

   for (; it.isValid (); ++it) {
      const SxString key_ = it->getProperty (".key").getString ();
      if (key_ == "")  res.prepend ("-");
      else             res.prepend (key_);
   }
   return res;
}

ssize_t SxSymbol::getSize () const
{
   SX_TRACE ();

   if (type != Type::List)
      SX_THROW ( "SxSymbol", "InvalidOperation",
                 SxString ("Expected type 'list', got '")
               + SxParserAst::getTypeStr (type) + "'");

   return dataG->begin (nodeId).getSizeOut ();
}

SxPtr<SxGraph<SxGProps> > SxSymbol::getGraphPtr () const
{
   SX_TRACE ();
   return dataG;
}

bool SxSymbol::isRoot () const
{
   SX_TRACE ();

   if (!isValid ())
      SX_THROW ("SxSymbol", "InvalidOperation",
                SxString ("Symbol is not valid"));

   return (dataG->begin(nodeId).getSizeIn () == 0);
}

void SxSymbol::write (const SxString &filename) const
{
   SX_TRACE ();
   try {
      SxFileIO f(filename, "w");
      SxSymbol (dataG, dataG->begin ().getIdx ()).print (f, "\n", 0, 3, true);
   } catch (SxException e) {
      SX_RETHROW (e, "SxSymbol", "IOError", SxString ("File write failed"));
   }
}

ssize_t SxSymbol::getIdx () const
{
   ssize_t res = SxUUIDv4 ().getHash ();

   while (dataG->containsNode (SxGProps (res)))  {
      res = SxUUIDv4 ().getHash ();
   }

   return res;
}

bool SxSymbol::isValidDocTxt (const SxString &docTxt)
{
   // regex to match SxDoc syntax
   SxRegex re("^\\/\\*[ \t]*SxDoc[ \t]*{[ \t]*[\r]?\n((.|\n|\r)*)}[ \t]*\\*\\/$");
   return (re.matchAll (docTxt).getSize () > 0);
}

} /* SxParserKit */
