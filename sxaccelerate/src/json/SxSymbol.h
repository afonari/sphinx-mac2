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

#ifndef _SX_SYMBOL_H_
#define _SX_SYMBOL_H_

#include <SxJSON.h>
#include <SxException.h>
#include <SxParserAst.h>
#include <SxGraph.h>
#include <SxGProps.h>
#include <SxGQuery.h>

/** \brief Wrapper for tree structured data

    \b SxSymbol = SPHInX Tree format Symbol

     SxSymbol class provides functionality
     to create, read, update or delete
     json data represented by SxGraph.

     SxSymbol class provides access
     functions to retrieve data out of
     file stored in tree structure that
     is represented by SxGraph.
 */

SX_EXCEPTION_CAT ("SxSymbol");
SX_EXCEPTION_TAG ("IOError");
SX_EXCEPTION_TAG ("InvalidInput");
SX_EXCEPTION_TAG ("InvalidOperation");

SX_EXCEPTION_TYPE ("SxSymbol", "IOError", "message", SxString);
SX_EXCEPTION_TYPE ("SxSymbol", "InvalidInput", "message", SxString);
SX_EXCEPTION_TYPE ("SxSymbol", "InvalidOperation", "message", SxString);

namespace SxParserKit {
class SX_EXPORT_JSON SxSymbol
{
   public:

      typedef typename SxGQExprBase::SelSet     SelSet;
      typedef typename SxParserAst::ElemType    Type;

      SxSymbol ();
      SxSymbol (const SxPtr<SxGraph<SxGProps> > &gPtr_, ssize_t nodeId_);
      SxSymbol (const SxSymbol &in) = delete;
      SxSymbol (SxSymbol &&in) = default;
     ~SxSymbol ();

      SxSymbol &operator= (const SxSymbol &in) = delete;
      SxSymbol &operator= (SxSymbol &&in) = default;

      // --- all to.. functions convert the current SxSymbol
      //     object to a specific type.
      SxString toString () const;
      int64_t  toInt () const;
      double   toDouble () const;
      bool     toBool () const;

      // --- functions to convert the SxSymbol
      //     object of type List to SxList of
      //     specific type.
      SxList<int64_t>   toIntList ()    const;
      SxList<SxString>  toStringList () const;
      SxList<double>    toDoubleList () const;
      SxList<SxSymbol>  toList ()       const;

      // --- check existence of an element
      bool hasElem (const SxString &key) const;

      // --- retrieve elem with given key
      //     inside the current SxSymbol object.
      SxSymbol getElem (const SxString &key, bool isRequired = true) const;
      // get an element at given index from current list
      SxSymbol getElem (ssize_t idx) const;

      // --- retrieve all elements of Group type
      SxList<SxSymbol> getGroupList (bool recursive = false) const;
      // --- retrieve all elements of Group type and given 'key_'
      SxList<SxSymbol> getGroupList (const SxString &key_,
                                     bool recursive = false) const;
      // --- retrieve all elements of List type
      SxList<SxSymbol> getArrayList (bool recursive = false) const;
      // --- retrieve all elements of List type and given 'key_'
      SxList<SxSymbol> getArrayList (const SxString &key_,
                                     bool recursive = false) const;

      SxString  getKey ()       const;
      SxVariant getValue ()     const;
      Type      getType ()      const;
      SxString  getTag  ()      const;
      SxString  getDocTxt ()    const;
      SxString  getDocTxtTag () const;

      // is current element valid
      bool isValid () const;

      // is current element the root object
      bool isRoot () const;
      // update an existing list element
      void setElem (ssize_t idx, const SxVariant &val);
      // update an existing group element
      void setElem (const SxString &key, const SxVariant &val);
      // set/change a type
      void setType (Type type);
      // set key of current element
      void setKey (const SxString &key);
      // set value of current element; might require type change
      void setValue (const SxVariant &val);
      // set SxDoc
      void setDocTxt (const SxString &str, const SxString &tag="");


      // get path to current element
      SxList<SxString> getPath () const;
      // get size of current list
      ssize_t getSize () const;
      // get underlying graph ptr
      SxPtr<SxGraph<SxGProps> > getGraphPtr () const;

      // finds the first element with given key in sub-tree
      SxSymbol find (const SxString &key) const;
      // find all elements with given key in sub-tree
      SxList<SxSymbol> findAll (const SxString &key) const;
      // get parent node except for root
      SxSymbol getParent () const;
      // append element of given type to array
      SxSymbol append (Type type_);
      // append element of given type to a group
      SxSymbol append (const SxString &key, Type type_);
      // append element to an array
      SxSymbol append (const SxVariant &val);
      // append given key/val to a group
      SxSymbol append (const SxString &key, const SxVariant &val);
      // remove given key from a group
      void remove (const SxString &key);
      // remove element from an array
      void remove (ssize_t idx);

      // write current json data to given file
      void write (const SxString &filename) const;
      // print json data to given stream
      void print (ostream &os=cout, const SxString &newLine="\n",
                  int lvl = 0, int tabs = 3,
                  bool isFileStream = false) const;

   protected:
      // check if valid doc text
      static bool isValidDocTxt (const SxString &docTxt);

      // --- retrieve all elements of given type
      SxList<SxSymbol> getTypedList (Type type_, bool recursive = false) const;
      SxList<SxSymbol> getTypedList (const SxString &key_, Type type_,
                                     bool recursive = false) const;

      // --- print in JSON format
      void printObject (ostream &os, const SxString &newLine,
                        int lvl, int tabs = 3,
                        bool isFileStream = false) const;
      void printArray (ostream &os, const SxString &newLine,
                       int lvl, int tabs = 3,
                       bool isFileStream = false) const;
      void printInt (ostream &os, const SxString &newLine,
                     int lvl, int tabs = 3,
                     bool isFileStream = false) const;
      void printDouble (ostream &os, const SxString &newLine,
                        int lvl, int tabs = 3,
                        bool isFileStream = false) const;
      void printString (ostream &os, const SxString &newLine,
                        int lvl, int tabs = 3,
                        bool isFileStream = false) const;
      void printBool (ostream &os, const SxString &newLine,
                      int lvl, int tabs = 3,
                      bool isFileStream = false) const;

      SxString name; // key
      Type type; // type of value
      SxPtr<SxGraph<SxGProps> > dataG;
      ssize_t nodeId; // graph node Id

      void removeRecursive (ssize_t nodeId_);
      ssize_t getIdx () const;
};

} /* SxParserKit */

#endif /* _SX_SYMBOL_H_ */
