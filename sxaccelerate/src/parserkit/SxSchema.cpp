#include <SxSchema.h>
#include <SxRegex.h>

SxSchema::SxSchema ()
   : schemaG (),
     dataG (),
     topLevelDefs (false),
     allowUndeclared (false)
{
   // empty
}

SxSchema::~SxSchema ()
{
   // empty
}

SxSchema::SxSchema (const SxPtr<SxGraph<SxGProps> > &schemaG_)
{
   SX_CHECK (schemaG_.getPtr ());

   schemaG  = schemaG_;
   simplifySchemaAst ();
}

bool SxSchema::getBool (const SxGraph<SxGProps>::Iterator &sIt,
                        const SxString &key,
                        bool default_)
{
   auto it = get (sIt, ".key", key);
   if (it.isValid ())
      return it->getProperty (".val").toBool ();
   else
      return default_;
}

int64_t SxSchema::getInt (const SxGraph<SxGProps>::Iterator &sIt,
                          const SxString &key,
                          const int64_t &default_)
{
   auto it = get (sIt, ".key", key);
   if (it.isValid ())
      return it->getProperty (".val").toInt ();
   else
      return default_;
}

double SxSchema::getDouble (const SxGraph<SxGProps>::Iterator &sIt,
                            const SxString &key,
                            const double &default_)
{
   auto it = get (sIt, ".key", key);
   if (it.isValid ())
      return it->getProperty (".val").toDouble ();
   else
      return default_;
}

SxString SxSchema::getString (const SxGraph<SxGProps>::Iterator &sIt,
                              const SxString &key,
                              const SxString &default_)
{
   auto it = get (sIt, ".key", key);
   if (it.isValid ())
      return it->getProperty (".val").toString ();
   else
      return default_;
}

void SxSchema::validationError (const SxString &msg,
                                const SxString &sTag,
                                const SxString &dTag)
{
   SxString errMsg = "";
   if (dTag != "")
      errMsg += dTag + ": ";
   errMsg += msg;
   if (sTag != "")
      errMsg += ", declared at " + sTag;
   SX_THROW ("SchemaValidation", "ValidationError",
             errMsg);
}

bool SxSchema::validateArray (const SxGraph<SxGProps>::Iterator &sIt,
                              const SxGraph<SxGProps>::Iterator &dIt)
{
   SX_CHECK (sIt.isValid ());

   int64_t n        = getInt (sIt, "nItems", max<int64_t>());
   int64_t maxItems = getInt (sIt, "maxItems", max<int64_t>());
   int64_t minItems = getInt (sIt, "minItems", min<int64_t>());

   bool isOptional = getBool (sIt, "optional", false);

   int64_t nItemsFound = dIt.getSizeOut ();

   if (nItemsFound == 0 && isOptional)  return true;

   if (n != max<int64_t>() && nItemsFound > n)  {
      validationError (  SxString("too many elements in list ")
                       + ", expected " + n
                       + " got " + nItemsFound,
                       sIt->getProperty (".val").getTag (),
                       dIt->getProperty (".val").getTag ());
   }
   if (maxItems != max<int64_t>() && nItemsFound > maxItems)  {
      validationError (  SxString("too many elements in list")
                       + ", expected less or equal to "
                       + maxItems + " got " + nItemsFound,
                       sIt->getProperty (".val").getTag (),
                       dIt->getProperty (".val").getTag ());
   }


   if (n != max<int64_t>() && nItemsFound < n)  {
      validationError (  SxString("too few elements in list, expected ")
                       + n + " got " + nItemsFound,
                       sIt->getProperty (".val").getTag (),
                       dIt->getProperty (".val").getTag ());
   }

   if (minItems != min<int64_t>() && nItemsFound < minItems)  {
      validationError (  SxString("too few elements in list, ")
                       + "expected greater or equal to "
                       + minItems + " got " + nItemsFound,
                       sIt->getProperty (".val").getTag (),
                       dIt->getProperty (".val").getTag ());
   }

   auto itemsIt = get (sIt, ".key", "items");

   Type type = (Type)getInt (itemsIt, "type", (int)Type::Undefined);

   switch (type)  {

      case Type::Group:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);

            if (elemIt->getProperty (".type").getInt () != Type::Group)  {
               Type t = (Type)elemIt->getProperty (".type").getInt ();
               SxString typeStr = SxParserAst::getTypeStr (t);
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'group', got '"
                                + typeStr + "'",
                                itemsIt->getProperty (".val").getTag (),
                                elemIt->getProperty (".val").getTag ());
            }

            validateObj (itemsIt, &elemIt);
         }
      }
      break;
      case Type::List:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);

            if (elemIt->getProperty (".type").getInt () != Type::List)  {
               Type t = (Type)elemIt->getProperty (".type").getInt ();
               SxString typeStr = SxParserAst::getTypeStr (t);
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'list', got '"
                                + typeStr + "'",
                                itemsIt->getProperty (".val").getTag (),
                                elemIt->getProperty (".val").getTag ());
            }

            validateArray (itemsIt, elemIt);
         }
      }
      break;
      case Type::Int:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);
            const SxVariant &val = elemIt->getProperty (".val");

            if (   val.getType () != SxVariantType::Int
                && val.getType () != SxVariantType::Double)  {
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'int', got '"
                                + val.getTypeName () + "'",
                                itemsIt->getProperty (".val").getTag (),
                                val.getTag ());
            }

            if (!itemsIt->getProperty (".val").matchLimits (val))  {
               int64_t minVal = 0;
               int64_t maxVal = 0;
               itemsIt->getProperty (".val").getLimits (&minVal, &maxVal);
               validationError (  SxString("int value ")
                                + val.toInt () + " is out of range ["
                                + minVal + ", "
                                + maxVal + "]",
                                itemsIt->getProperty (".val").getTag (),
                                val.getTag ());
            }

         }
      }
      break;
      case Type::Double:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);
            const SxVariant &val = elemIt->getProperty (".val");

            if (   val.getType () != SxVariantType::Int
                && val.getType () != SxVariantType::Double)
            {
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'real', got '"
                                + val.getTypeName () + "'",
                                itemsIt->getProperty (".val").getTag (),
                                val.getTag ());
            }

            if (!itemsIt->getProperty (".val").matchLimits (val))  {
               double minVal = 0;
               double maxVal = 0;
               itemsIt->getProperty (".val").getLimits (&minVal, &maxVal);
               validationError (  SxString("real value ")
                                + val.toDouble () + " is out of range ["
                                + minVal + ", "
                                + maxVal + "]",
                                itemsIt->getProperty (".val").getTag (),
                                val.getTag ());
            }
         }
      }
      break;
      case Type::Bool:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);
            SxVariant &val = elemIt->getProperty (".val");

            // --- correction for symboltable issue
            if (val.getType () == SxVariantType::Int)  {
               int64_t v = elemIt->getProperty (".val").toInt ();
               if (v > 0)  elemIt->setProperty (".val", true);
               else        elemIt->setProperty (".val", false);
               elemIt->setProperty (".type", (int)Type::Bool);
            }

            if (val.getType () != SxVariantType::Bool)  {
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'bool', got '"
                                + val.getTypeName () + "'",
                                itemsIt->getProperty (".val").getTag (),
                                val.getTag ());
            }

         }
      }
      break;
      case Type::String:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);
            const SxVariant &val = elemIt->getProperty (".val");

            if (val.getType () != SxVariantType::String)
            {
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'string', got '"
                                + val.getTypeName () + "'",
                                itemsIt->getProperty (".val").getTag (),
                                val.getTag ());
            }

            if (!itemsIt->getProperty (".val").matchLimits (val))  {
               ssize_t minLen = 0;
               ssize_t maxLen = 0;
               itemsIt->getProperty (".val").getLimits (&minLen, &maxLen);
               ssize_t strLen = val.toString ().getSize ();
               if (strLen < minLen || strLen > maxLen)  {
                  validationError (  SxString("string of size ")
                                   + val.toString ().getSize ()
                                   + " is out of range ["
                                   + minLen + ", "
                                   + maxLen + "]",
                                   itemsIt->getProperty (".val").getTag (),
                                   val.getTag ());
               }  else  {
                  SxString regexHint = getString (itemsIt, "regexHint", "");
                  validationError (  SxString ("unexpected format for list element")
                                   + ", expected '"
                                   + regexHint + "'",
                                   itemsIt->getProperty (".val").getTag (),
                                   val.getTag ());
               }
            }
         }
      }
      break;
      case Type::Enum:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);
            if (   elemIt->getProperty (".val").getType ()
                != SxVariantType::String)
            {
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'string(enum)', got '"
                                + elemIt->getProperty (".val").getTypeName ()
                                + "'",
                                itemsIt->getProperty (".val").getTag (),
                                elemIt->getProperty (".val").getTag ());
            }

            SxString val = elemIt->getProperty (".val").getString ();

            SxList<SxString> options = getString (itemsIt, "val", "").tokenize(",");
            bool found = false;
            for (auto oIt = options.begin (); oIt.isValid (); ++oIt)  {
               if (*oIt == val)  {
                  found = true;
                  break;
               }
            }

            if (found == false)  {
               validationError (  SxString ("value for list element")
                                + " not part of Enum values",
                                itemsIt->getProperty (".val").getTag (),
                                elemIt->getProperty (".val").getTag ());
            }
         }
      }
      break;
      case Type::Vector:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);
            if (elemIt->getProperty (".type").getInt () != Type::List)  {
               Type t = (Type)elemIt->getProperty (".type").getInt ();
               SxString typeStr = SxParserAst::getTypeStr (t);
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'list(vector)', got '"
                                + typeStr + "'",
                                itemsIt->getProperty (".val").getTag (),
                                elemIt->getProperty (".val").getTag ());
            }

            validateVector (itemsIt, elemIt);
         }
      }
      break;
      case Type::Matrix:  {
         for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
            auto elemIt = dIt.out (i);
            if (elemIt->getProperty (".type").getInt () != Type::List)  {
               Type t = (Type)elemIt->getProperty (".type").getInt ();
               SxString typeStr = SxParserAst::getTypeStr (t);
               validationError (  SxString("type mismatch for list element")
                                + ", expected 'list(matrix)', got '"
                                + typeStr + "'",
                                itemsIt->getProperty (".val").getTag (),
                                elemIt->getProperty (".val").getTag ());
            }

            validateMatrix (itemsIt, elemIt);
         }
      }
      break;
      default:  {
         validationError ("list element type is not recognized",
                          "",
                          itemsIt->getProperty (".val").getTag ());
      }
   }

   return true;
}

bool SxSchema::validateVector (const SxGraph<SxGProps>::Iterator &sIt,
                               const SxGraph<SxGProps>::Iterator &dIt)
{
   SX_CHECK (sIt.isValid ());

   ssize_t nElems = dIt.getSizeOut ();

   int64_t dim = getInt (sIt, "dim", max<int64_t>());

   if (dim != max<int64_t>() && dim != static_cast<int64_t>(nElems))  {
      validationError (  SxString("dimensions mismatch for vector, ")
                       + "must be equal to " + dim,
                       sIt->getProperty (".val").getTag (),
                       dIt->getProperty (".val").getTag ());
   }

   int64_t minDim = getInt (sIt, "minDim", min<int64_t>());

   if (minDim != min<int64_t>() && static_cast<int64_t>(nElems) < minDim)  {
      validationError (  SxString("dimensions mismatch for vector, ")
                       + "must be greater or equal to "
                       + minDim + ", got " + nElems,
                       sIt->getProperty (".val").getTag (),
                       dIt->getProperty (".val").getTag ());
   }

   int64_t maxDim = getInt (sIt, "maxDim", max<int64_t>());

   if (maxDim != max<int64_t>() && static_cast<int64_t>(nElems) > maxDim)  {
      validationError (  SxString("dimensions mismatch for vector, ")
                       + "must be less or equal to "
                       + maxDim + ", got " + nElems,
                       sIt->getProperty (".val").getTag (),
                       dIt->getProperty (".val").getTag ());
   }

   return true;
}

bool SxSchema::validateMatrix (const SxGraph<SxGProps>::Iterator &sIt,
                               const SxGraph<SxGProps>::Iterator &dIt)
{
   SX_CHECK (sIt.isValid ());


   int64_t dimX = max<int64_t>(), dimY = max<int64_t>();

   auto dimsIt = get (sIt, ".key", "dims");
   SX_CHECK (dimsIt.isValid ());
   SX_CHECK (dimsIt->getProperty (".type").getInt () == Type::List);

   if (   dimsIt.isValid ()
       && dimsIt->getProperty (".type").getInt () == Type::List)  {

      if (dimsIt.getSizeOut () > 0)
         dimX = dimsIt.out (0)->getProperty (".val").toInt ();

      if (dimsIt.getSizeOut () > 1)
         dimY = dimsIt.out (1)->getProperty (".val").toInt ();

      SX_CHECK (dimX != 0 && dimY != 0);

      ssize_t nElems = dIt.getSizeOut ();

      if (dimX != max<int64_t>() && dimX != static_cast<int64_t>(nElems))  {
         validationError (  SxString("number of rows ")
                          + "of matrix "
                          + "must be equal to " + dimX
                          + ", got " + nElems,
                          sIt->getProperty (".val").getTag (),
                          dIt->getProperty (".val").getTag ());
      }

      nElems = 0;

      for (ssize_t i = 0; i < dIt.getSizeOut (); ++i)  {
         auto eIt = dIt.out (i);
         SX_CHECK (eIt->getProperty (".type").toInt () == Type::List);
         nElems = eIt.getSizeOut ();
         if (dimY != static_cast<int64_t>(nElems))  {
            validationError (  SxString("number of columns ")
                             + " of matrix "
                             + "must be equal to " + dimY
                             + ", got " + nElems,
                             sIt->getProperty (".val").getTag (),
                             eIt->getProperty (".val").getTag ());
         }
      }
   }
   return true;
}

bool SxSchema::validateObj (const SxGraph<SxGProps>::Iterator &sIt,
                            SxGraph<SxGProps>::Iterator *dIt)
{
   SX_CHECK (dIt);
   SX_CHECK (sIt.isValid ());

   auto dataIt = *dIt;

   // --- any unexpected items
   if (!allowUndeclared)  {
      for (ssize_t m = 0; m < dataIt.getSizeOut (); ++m)  {
         bool found = false;
         SxString key = dataIt.out (m)->getProperty (".key").getString ();
         for (ssize_t n = 0; n < sIt.getSizeOut (); ++n)  {
            SxString k = sIt.out (n)->getProperty (".key").getString ();
            if (key == k)  {
               found = true;
               break;
            }
         }
         if (!found)  {
            validationError (  "undeclared key '"
                             + key + "' in group '"
                             + sIt->getProperty (".key").getString ()
                             + "'",
                             sIt->getProperty (".val").getTag (),
                             dataIt.out (m)->getProperty (".key").getTag ());
         }
      }
   }

   SxMap<SxString, ssize_t> nameToIdx;
   SxArray<uint8_t> isFound;
   isFound.resize (sIt.getSizeOut ());
   if (isFound.getSize () > 0)  isFound.set (0);
   for (ssize_t i = 0; i < sIt.getSizeOut (); ++i)  {
      auto chIt = sIt.out (i);

      if (chIt->getProperty (".type").getInt () != Type::Group)
         continue;

      Type type = (Type)getInt (chIt, "type", (int)Type::Undefined);

      const SxString &key = chIt->getProperty (".key").getString ();

      nameToIdx (key) = i;

      switch (type)  {

         case Type::Group:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  if (dchIt->getProperty (".type").getInt () != Type::Group)  {
                     Type t = (Type)dchIt->getProperty (".type").getInt ();
                     SxString typeStr = SxParserAst::getTypeStr (t);
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'group', got '"
                                      + typeStr + "'",
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".val").getTag ());
                  }
                  prevId = k;
                  isFound(i) = 1;
                  validateObj (chIt, &dchIt);
               }
            }
         }
         break;
         case Type::List:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  if (dchIt->getProperty (".type").getInt () != Type::List)  {
                     Type t = (Type)dchIt->getProperty (".type").getInt ();
                     SxString typeStr = SxParserAst::getTypeStr (t);
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'list', got '"
                                      + typeStr + "'",
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".val").getTag ());
                  }
                  prevId = k;
                  isFound(i) = 1;
                  validateArray (chIt, dchIt);
               }
            }

         }
         break;
         case Type::Vector:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  if (dchIt->getProperty (".type").getInt () != Type::List)  {
                     Type t = (Type)dchIt->getProperty (".type").getInt ();
                     SxString typeStr = SxParserAst::getTypeStr (t);
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'list(vector)', got '"
                                      + typeStr + "'",
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".val").getTag ());
                  }
                  prevId = k;
                  isFound(i) = 1;
                  validateVector (chIt, dchIt);
               }
            }

         }
         break;
         case Type::Matrix:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  if (dchIt->getProperty (".type").getInt () != Type::List)  {
                     Type t = (Type)dchIt->getProperty (".type").getInt ();
                     SxString typeStr = SxParserAst::getTypeStr (t);
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'list(matrix)', got '"
                                      + typeStr + "'",
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".val").getTag ());
                  }
                  prevId = k;
                  isFound(i) = 1;
                  validateMatrix (chIt, dchIt);
               }
            }

         }
         break;
         case Type::Enum:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  isFound(i) = 1;
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  if (   dchIt->getProperty (".val").getType ()
                      != SxVariantType::String)
                  {
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'string(enum)', got '"
                                      + dchIt->getProperty (".val").getTypeName ()
                                      + "'",
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".val").getTag ());
                  }

                  SxString val = dchIt->getProperty (".val").getString ();

                  SxList<SxString> options = getString (chIt, "val", "").tokenize(",");
                  bool found = false;
                  for (auto oIt = options.begin (); oIt.isValid (); ++oIt)  {
                     if (*oIt == val)  {
                        found = true;
                        break;
                     }
                  }

                  if (found == false)  {
                     validationError (  SxString ("value for key '" + key
                                      + "' not part of Enum values"),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".val").getTag ());
                  }
                  prevId = k;
               }
            }

         }
         break;
         case Type::Int:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  isFound(i) = 1;
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  const SxVariant &val = dchIt->getProperty (".val");

                  if (   val.getType () != SxVariantType::Int
                      && val.getType () != SxVariantType::Double)  {
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'int', got '"
                                      + val.getTypeName () + "'",
                                      chIt->getProperty (".val").getTag (),
                                      val.getTag ());
                  }

                  if (!chIt->getProperty (".val").matchLimits (val))  {
                     int64_t minVal = 0;
                     int64_t maxVal = 0;
                     chIt->getProperty (".val").getLimits (&minVal, &maxVal);
                     validationError (  SxString("int value ")
                                      + val.toInt () + " is out of range ["
                                      + minVal + ", "
                                      + maxVal + "]",
                                      chIt->getProperty (".val").getTag (),
                                      val.getTag ());
                  }
                  prevId = k;
               }
            }

         }
         break;
         case Type::Double:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  isFound(i) = 1;
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  const SxVariant &val = dchIt->getProperty (".val");

                  if (   val.getType () != SxVariantType::Int
                      && val.getType () != SxVariantType::Double)
                  {
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'real', got '"
                                      + val.getTypeName () + "'",
                                      chIt->getProperty (".val").getTag (),
                                      val.getTag ());
                  }

                  if (!chIt->getProperty (".val").matchLimits (val))  {
                     double minVal = 0;
                     double maxVal = 0;
                     chIt->getProperty (".val").getLimits (&minVal, &maxVal);
                     validationError (  SxString("real value ")
                                      + val.toDouble () + " is out of range ["
                                      + minVal + ", "
                                      + maxVal + "]",
                                      chIt->getProperty (".val").getTag (),
                                      val.getTag ());
                  }
                  prevId = k;
               }
            }

         }
         break;
         case Type::Bool:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  isFound(i) = 1;
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  SxVariant &val = dchIt->getProperty (".val");

                  // --- correction for symboltable issue
                  if (val.getType () == SxVariantType::Int)  {
                     int64_t v = dchIt->getProperty (".val").toInt ();
                     if (v > 0)  dchIt->setProperty (".val", true);
                     else        dchIt->setProperty (".val", false);
                     dchIt->setProperty (".type", (int)Type::Bool);
                  }

                  if (val.getType () != SxVariantType::Bool)  {
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'bool', got '"
                                      + val.getTypeName () + "'",
                                      chIt->getProperty (".val").getTag (),
                                      val.getTag ());
                  }
                  prevId = k;
               }
            }

         }
         break;
         case Type::String:  {
            ssize_t prevId = -1;
            for (ssize_t k = 0; k < dataIt.getSizeOut (); ++k)  {
               auto dchIt = dataIt.out (k);
               if (dchIt->getProperty (".key").getString () == key)  {
                  isFound(i) = 1;
                  if (prevId != -1)  {
                     validationError (  "duplicate key '"
                                      + key
                                      + "' previously defined at "
                                      + dataIt.out (prevId)->getProperty (".key").getTag (),
                                      chIt->getProperty (".val").getTag (),
                                      dchIt->getProperty (".key").getTag ());
                  }

                  const SxVariant &val = dchIt->getProperty (".val");

                  if (val.getType () != SxVariantType::String)
                  {
                     validationError (  "type mismatch for key '" + key
                                      + "', expected 'string', got '"
                                      + val.getTypeName () + "'",
                                      chIt->getProperty (".val").getTag (),
                                      val.getTag ());
                  }

                  if (!chIt->getProperty (".val").matchLimits (val))  {
                     ssize_t minLen = 0;
                     ssize_t maxLen = 0;
                     chIt->getProperty (".val").getLimits (&minLen, &maxLen);
                     ssize_t strLen = val.toString ().getSize ();
                     if (strLen < minLen || strLen > maxLen)  {
                        validationError (  SxString("string of size ")
                                         + val.toString ().getSize ()
                                         + " is out of range ["
                                         + minLen + ", "
                                         + maxLen + "]",
                                         chIt->getProperty (".val").getTag (),
                                         val.getTag ());
                     }  else  {
                        SxString regexHint = getString (chIt, "regexHint", "");
                        validationError (  SxString ("unexpected format for key '")
                                         + key + "', expected '"
                                         + regexHint + "'",
                                         chIt->getProperty (".val").getTag (),
                                         val.getTag ());
                     }
                  }
                  prevId = k;
               }
            }

         }
         break;
         default:  {
            validationError ("type is not recognized",
                             "",
                             chIt->getProperty (".val").getTag ());
         }

      }

      bool opt = getBool (chIt, "optional", false);

      // --- if not found and is optional, set to 2
      if (isFound(i) == 0 && opt == true)  {
         isFound(i) = 2;
      }

      // ----------------------------------------------------

      if (isFound(i) == 1)  {
         // --- check for needed items
         SxString needs = getString (chIt, "needs", "");

         if (needs != "")  {
            SxList<SxString> needList = needs.tokenize (",");
            bool found = true;
            for (auto nIt = needList.begin (); nIt.isValid (); ++nIt)
            {
               if (!(get(dataIt,".key", *nIt).isValid()))  {
                  found = false;
                  break;
               }
            }

            if (!found)  {
               validationError (  "Presence of item(s) '"
                                + needs
                                + "' is required by '"
                                + chIt->getProperty (".key").getString ()
                                + "'",
                                chIt->getProperty (".val").getTag (),
                                dataIt->getProperty (".val").getTag ()
                               );
            }
         }

         // --- check for excluded items
         SxString nots = getString (chIt, "not", "");

         if (nots != "")  {
            SxList<SxString> notList = nots.tokenize (",");
            bool found = false;
            for (auto nIt = notList.begin (); nIt.isValid (); ++nIt)
            {
               if ((get(dataIt,".key", *nIt).isValid()))  {
                  found = true;
                  break;
               }
            }

            if (found)  {
               validationError (  "Presence of item(s) '"
                                + nots
                                + "' is excluded by '"
                                + chIt->getProperty (".key").getString ()
                                + "'",
                                chIt->getProperty (".val").getTag (),
                                dataIt->getProperty (".val").getTag ()
                               );
            }
         }
      }

   }

   for (ssize_t i = 0; i < sIt.getSizeOut (); ++i)  {
      auto chIt = sIt.out (i);

      if (chIt->getProperty (".type").getInt () != Type::Group)
         continue;

      const SxString &key = chIt->getProperty (".key").getString ();
      SxString xorStr = getString (chIt, "xor", "");

      if (xorStr != "")  {
         SxList<SxString> xorList = xorStr.tokenize (',');
         bool xorFound = false;
         for (auto xorIt = xorList.begin (); xorIt != xorList.end (); ++xorIt)
         {
            if (   isFound(nameToIdx (*xorIt)) == 1
                && *xorIt != key)  {
               xorFound = true;
            }
         }

         if (xorFound == false && isFound (i) == 0)  {
            validationError (  "Missing one of '"
                             + xorStr + "'",
                             chIt->getProperty (".val").getTag (),
                             dataIt->getProperty (".val").getTag ());
         }
         // --- if curr key is found and one of xor's is also found
         if (xorFound == true && isFound (i) == 1)  {
            validationError (xorStr + "' are mutually excluded",
                             chIt->getProperty (".val").getTag (),
                             dataIt->getProperty (".val").getTag ());
         }
      }  else  {
         // --- if not optional
         if (isFound (i) == 0)  {
            // --- item not found
            validationError (  "key '" + key
                             + "' not found in '"
                             + dataIt->getProperty (".key").getString ()
                             + "'", sIt->getProperty (".val").getTag (),
                             dataIt->getProperty (".val").getTag ());
         }
      }
   }

   return true;
}

bool SxSchema::validate (const SxPtr<SxGraph<SxGProps> > &dataG_)
{
   SX_CHECK (dataG_.getPtr ());

   dataG = dataG_;

   auto dIt = dataG->begin ();
   auto sIt = schemaG->begin ();

   Type t = (Type)getInt (sIt, "type", Type::Undefined);
   if (t == Type::Group)  validateObj (sIt, &dIt);
   else                   validateArray (sIt, dIt);

   return true;
}

SxGraph<SxGProps>::Iterator
SxSchema::get (const SxGraph<SxGProps>::Iterator &it,
               const SxString &key,
               const SxVariant &val)
{
   SX_CHECK (key != "");
   SX_CHECK (it.isValid ());

   for (ssize_t i = 0; i < it.getSizeOut (); ++i)
   {
      auto chIt = it.out (i);
      if (chIt->hasProperty (key))  {
         if (!val.isInitialized () || chIt->getProperty (key) == val)  {
            return chIt;
         }
      }
   }
   return SxGraph<SxGProps>::Iterator ();
}

void SxSchema::simplifyInt (SxGraph<SxGProps>::Iterator &sIt)
{
   int64_t minVal = getInt (sIt, "min", min<int64_t>());
   int64_t maxVal = getInt (sIt, "max", max<int64_t>());

   sIt->getProperty (".val").setType (SxVariantType::Int);
   sIt->getProperty (".val").set (minVal);
   sIt->getProperty (".val").setLimits (minVal, maxVal);
   SX_CHECK (sIt->getProperty (".val").getType() == SxVariantType::Int);
}

void SxSchema::simplifyDouble (SxGraph<SxGProps>::Iterator &sIt)
{
   double minVal = getDouble (sIt, "min", min<double>());
   double maxVal = getDouble (sIt, "max", max<double>());

   sIt->getProperty (".val").setType (SxVariantType::Double);
   sIt->getProperty (".val").set (minVal);
   sIt->getProperty (".val").setLimits (minVal, maxVal);
   SX_CHECK (sIt->getProperty (".val").getType() == SxVariantType::Double);
}

void SxSchema::simplifyString (SxGraph<SxGProps>::Iterator &sIt)
{
   ssize_t minVal = (ssize_t)getInt (sIt, "min", min<ssize_t>());
   ssize_t maxVal = (ssize_t)getInt (sIt, "max", max<ssize_t>());

   sIt->getProperty (".val").setType (SxVariantType::String);
   sIt->getProperty (".val").set (SxString (""));
   sIt->getProperty (".val").setLimits (minVal, maxVal);

   SxString regex = getString (sIt, "regex", "");
   if (regex != "")  {
      try  {
         sIt->getProperty (".val").setRegex (regex);
      }  catch (SxException e)  {
         auto regexIt = get (sIt, ".key", "regex");
         validationError (e.toString (), "",
                          regexIt->getProperty (".val").getTag ());
      }
   }

   SX_CHECK (sIt->getProperty (".val").getType() == SxVariantType::String);
}

void SxSchema::simplifyArray (SxGraph<SxGProps>::Iterator &sIt)
{
   auto itemsIt = get (sIt, ".key", "items");
   SX_CHECK (itemsIt.isValid ());

   auto typeIt = get (itemsIt, ".key", "type");
   SxString typeStr = typeIt->getProperty (".val").toString ();
   Type itemType = (Type)SxParserAst::getTypeId (typeStr);
   typeIt->setProperty (".val", (int)itemType);

   if (itemType == Type::Group)  {
      simplifyObj (itemsIt);
   }  else if (itemType == Type::List)  {
      simplifyArray (itemsIt);
   }  else if (itemType == Type::Int)  {
      simplifyInt (itemsIt);
   }  else if (itemType == Type::Double)  {
      simplifyDouble (itemsIt);
   }  else if (itemType == Type::String)  {
      simplifyString (itemsIt);
   }

}

void SxSchema::simplifyObj (SxGraph<SxGProps>::Iterator &sIt)
{
   for (ssize_t i = 0; i < sIt.getSizeOut (); ++i)  {
      auto chIt = sIt.out (i);

      if (chIt->getProperty (".type").getInt () == Type::Group)  {

         auto typeIt = get (chIt, ".key", "type");

         SxString typeStr = typeIt->getProperty (".val").toString ();

         Type elemType = (Type)SxParserAst::getTypeId (typeStr);
         typeIt->setProperty (".val", (int)elemType);

         if (elemType == Type::Group)  {
            simplifyObj (chIt);
         }  else if (elemType == Type::Int)  {
            simplifyInt (chIt);
         }  else if (elemType == Type::Double)  {
            simplifyDouble (chIt);
         }  else if (elemType == Type::String)  {
            simplifyString (chIt);
         }  else if (elemType == Type::List)  {
            simplifyArray (chIt);
         }
      }

   }
}

void SxSchema::simplifySchemaAst ()
{
   SX_CHECK (schemaG.getPtr ());

   auto sIt = schemaG->begin ();

   topLevelDefs = getBool (sIt, "topLevelDefs", false);

   allowUndeclared = getBool (sIt, "allowUndeclared", false);

   auto typeIt = get (sIt, ".key", "type");

   Type t = (Type)SxParserAst::getTypeId (typeIt->getProperty (".val").toString ());
   typeIt->setProperty (".val", (int)t);

   if (t == Type::Group)  simplifyObj (sIt);
}

