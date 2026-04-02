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

#include <SxVariant.h>
#include <SxException.h>
#include <SxMacroLib.h>
#include <SxSerializer.h>

void SxVariant::setType (int type_)
{
   SX_TRACE ();
   if (type != type_)  {
      reset ();
      type = type_;
      switch (type)  {
         case SxVariantType::Undefined:
         case SxVariantType::Any:
            dataType = SxVariantType::Undefined;
            break;
         case SxVariantType::Int:
            dataType = SxVariantType::Int;
            data   = (void *)new int64_t;
            minPtr = (void *)new int64_t;
            maxPtr = (void *)new int64_t;
            *(int64_t *)data   = 0;
            *(int64_t *)minPtr = 0;
            *(int64_t *)maxPtr = 0;
            break;
         case SxVariantType::Double:
            dataType = SxVariantType::Double;
            data   = (void *)new double;
            minPtr = (void *)new double;
            maxPtr = (void *)new double;
            *(double *)data   = 0.0;
            *(double *)minPtr = 0.0;
            *(double *)maxPtr = 0.0;
            break;
         case SxVariantType::Bool:
            dataType = SxVariantType::Bool;
            data   = (void *)new bool;
            minPtr = (void *)new uint8_t;
            maxPtr = (void *)new uint8_t;
            *(bool *)data      = false;
            *(uint8_t *)minPtr = 0;
            *(uint8_t *)maxPtr = 0;
            break;
         case SxVariantType::String:
            dataType = SxVariantType::String;
            data   = (void *)new SxString;
            minPtr = (void *)new char;
            maxPtr = (void *)new char;
            *(char *)minPtr = 0;
            *(char *)maxPtr = 0;
            break;
         case SxVariantType::List:
            dataType = SxVariantType::List;
            data   = (void *)new SxList<SxVariant>;
            minPtr = (void *)new int;
            maxPtr = (void *)new int;
            *(int *)minPtr = SxVariantType::Undefined;
            *(int *)maxPtr = SxVariantType::Undefined;
            initialized = true;
            break;
         case SxVariantType::StringArray:
            dataType = SxVariantType::StringArray;
            data = (void *)new SxArray<SxString>;
            initialized = true;
            break;
         case SxVariantType::IntArray:
            dataType = SxVariantType::IntArray;
            data = (void *)new SxArray<int64_t>;
            initialized = true;
            break;
         case SxVariantType::FloatArray:
            dataType = SxVariantType::FloatArray;
            data = (void *)new SxArray<float>;
            initialized = true;
            break;
         case SxVariantType::DoubleArray:
            dataType = SxVariantType::DoubleArray;
            data = (void *)new SxArray<double>;
            initialized = true;
            break;
         default: SX_EXIT;
      }
      resetLimits ();
   }
}

void SxVariant::reset ()
{
   SX_TRACE ();
   if (data)  {
      switch (dataType)  {
         case SxVariantType::Undefined:
         case SxVariantType::Any:
            SX_EXIT;
            break;
         case SxVariantType::Int:
            delete ((int64_t *)data);
            delete ((int64_t *)minPtr);
            delete ((int64_t *)maxPtr);
            break;
         case SxVariantType::Double:
            delete ((double *)data);
            delete ((double *)minPtr);
            delete ((double *)maxPtr);
            break;
         case SxVariantType::Bool:
            delete ((bool *)data);
            delete ((uint8_t *)minPtr);
            delete ((uint8_t *)maxPtr);
            break;
         case SxVariantType::String:
            delete ((SxString *)data);
            delete ((char *)minPtr);
            delete ((char *)maxPtr);
            break;
         case SxVariantType::List:
            delete ((SxList<SxVariant> *)data);
            delete ((int *)minPtr);
            delete ((int *)maxPtr);
            break;
         case SxVariantType::StringArray:
            delete ((SxArray<SxString> *)data);
            break;
         case SxVariantType::IntArray:
            delete ((SxArray<int64_t> *)data);
            break;
         case SxVariantType::FloatArray:
            delete ((SxArray<float> *)data);
            break;
         case SxVariantType::DoubleArray:
            delete ((SxArray<double> *)data);
            break;
         default: SX_EXIT;
      }
      data = NULL;
      minPtr = NULL;
      maxPtr = NULL;
      initialized = false;
      type = SxVariantType::Undefined;
      dataType = SxVariantType::Undefined;
   }
}

void SxVariant::resetLimits ()
{
   SX_TRACE ();
   if (data)  {
      switch (dataType)  {
         case SxVariantType::Undefined:
         case SxVariantType::Any:
            break;
         case SxVariantType::Int:
            *(int64_t *)minPtr = SxVariantType::ScalarType<int>::minVal ();
            *(int64_t *)maxPtr = SxVariantType::ScalarType<int>::maxVal ();
            break;
         case SxVariantType::Double:
            *(double *)minPtr = SxVariantType::ScalarType<double>::minVal ();
            *(double *)maxPtr = SxVariantType::ScalarType<double>::maxVal ();
            break;
         case SxVariantType::Bool:
            *(uint8_t *)minPtr = SxVariantType::ScalarType<bool>::minVal ();
            *(uint8_t *)maxPtr = SxVariantType::ScalarType<bool>::maxVal ();
            break;
         case SxVariantType::String:
            *(char *)minPtr = SxVariantType::ScalarType<SxString>::minVal ();
            *(char *)maxPtr = SxVariantType::ScalarType<SxString>::maxVal ();
            break;
         case SxVariantType::List:
            *(int *)minPtr = SxVariantType::Undefined;
            *(int *)maxPtr = SxVariantType::Undefined;
            break;
         case SxVariantType::StringArray:
         case SxVariantType::IntArray:
         case SxVariantType::FloatArray:
         case SxVariantType::DoubleArray:
            break;
         default: SX_EXIT;
      }
      minElem = 0;
      maxElem = 0;
   }
}

void SxVariant::operator= (const SxVariant &in)
{
   SX_TRACE ();
   if (&in == this)  return;

   setType (in.type);

   switch (dataType) {
      case SxVariantType::Undefined:
      case SxVariantType::Any:
         break;
      case SxVariantType::Int:
         *(int64_t *)data   = *(int64_t *)in.data;
         *(int64_t *)minPtr = *(int64_t *)in.minPtr;
         *(int64_t *)maxPtr = *(int64_t *)in.maxPtr;
         break;
      case SxVariantType::Double:
         *(double *)data   = *(double *)in.data;
         *(double *)minPtr = *(double *)in.minPtr;
         *(double *)maxPtr = *(double *)in.maxPtr;
         break;
      case SxVariantType::Bool:
         *(bool *)data      = *(bool *)in.data;
         *(uint8_t *)minPtr = *(uint8_t *)in.minPtr;
         *(uint8_t *)maxPtr = *(uint8_t *)in.maxPtr;
         break;
      case SxVariantType::String:
         *(SxString *)data  = *(SxString *)in.data;
         *(char *)minPtr = *(char *)in.minPtr;
         *(char *)maxPtr = *(char *)in.maxPtr;
         break;
      case SxVariantType::List:
         *(SxList<SxVariant> *)data = *(SxList<SxVariant> *)in.data;
         *(int *)minPtr = *(int *)in.minPtr;
         *(int *)maxPtr = *(int *)in.maxPtr;
         break;
      case SxVariantType::StringArray:
         *(SxArray<SxString> *)data = in.toStringArray ();
         break;
      case SxVariantType::IntArray:
         *(SxArray<int64_t> *)data = in.toIntArray ();
         break;
      case SxVariantType::FloatArray:
         *(SxArray<float> *)data = in.toFloatArray ();
         break;
      case SxVariantType::DoubleArray:
         *(SxArray<double> *)data = in.toDoubleArray ();
         break;
      default: SX_EXIT;
   }

   initialized = in.initialized;
   minElem = in.minElem;
   maxElem = in.maxElem;
   regex   = in.regex;
   regexPattern = in.regexPattern;
   regexFlags   = in.regexFlags;

   tag = in.tag;
}

bool SxVariant::isValidString (const SxString &in) const
{
   SX_TRACE ();
   if (minElem && in.getSize () < minElem)
      return false;
   else if (maxElem && in.getSize () > maxElem)
      return false;
   if (regex.getPtr() != NULL)  {
      try  {
         return (regex->match(in).getSize() == 1);
      }  catch (SxException e)  {
         return false;
      }
   }

   return true;
}

void SxVariant::validate (const SxString &in) const
{
   SX_TRACE ();
   if (minElem && in.getSize () < minElem)
      SX_THROW("String is below minimal length " + SxString(minElem));
   else if (maxElem && in.getSize () > maxElem)
      SX_THROW("String exceeds maximum length " + SxString(maxElem));
}

bool SxVariant::hasStringLimits () const
{
   SX_TRACE ();
   if (regex.getPtr() != NULL || minElem || maxElem) return true;
   else                                              return false;
}

bool SxVariant::matchLimits (const SxVariant &in) const
{
   SX_TRACE ();
   using namespace SxVariantType;

   int typeB = in.getType ();

   if (type == Int && typeB == Int)  {
      int64_t val = in.getInt ();
      return !(val < minInt() || val > maxInt());
   }  else if (type == Int && typeB == Double)  {
      int64_t val = static_cast<int64_t>(in.getDouble ());
      return !(val < minInt() || val > maxInt());
   }  else if (type == Double && typeB == Double)  {
      double val = in.getDouble ();
      return !(val <  minDouble() || val > maxDouble());
   }  else if (type == Double && typeB == Int)  {
      double val = static_cast<double>(in.getInt ());
      return !(val <  minDouble() || val > maxDouble());
   }  else if (type == String && typeB == String)  {
      if (hasStringLimits ())
         return isValidString (in.toString ());
      else
         return true;
   }

   return false;
}

bool SxVariant::match (const SxVariant &in, bool withLimits_) const
{
   SX_TRACE ();
   using namespace SxVariantType;

   if (type == SxVariantType::Any)  {
      return true;
   } else if (type == in.getType ())  {
      if (withLimits_)  return matchLimits(in);
      else              return true;

   } else if (   (type == Int && in.getType () == Double)
              || (type == Double && in.getType () == Int))
   {
      if (withLimits_)  return matchLimits(in);
      else              return true;
   }

   return false;
}

SxVariant::CastType SxVariant::castDiff (const SxVariant &in) const
{
   SX_TRACE ();
   using namespace SxVariantType;
   if (type == in.type)       return  ExactCast;
   if (type == Any)           return  ExactCast;
   if (in.type == Undefined)  return  IncompatibleTypes;

   if (type == Int)  {
      switch (in.type)  {
         case Double: return ImplicitCast;
         case Bool:   return ImplicitCast;
         case String: return ExplicitCastRequired;
         case List:   return ExplicitCastRequired;
      }
   } else if (type == Double)  {
      switch (in.type)  {
         case Int:    return ExplicitCastRequired;
         case Bool:   return ImplicitCast;
         case String: return ExplicitCastRequired;
         case List:   return ExplicitCastRequired;
      }
   } else if (type == Bool)  {
      switch (in.type)  {
         case Int:    return ImplicitCast;
         case Double: return ImplicitCast;
         case String: return ExplicitCastRequired;
         case List:   return ExplicitCastRequired;
      }
   } else if (type == String)  {
      switch (in.type)  {
         case Int:    return  ExplicitCastRequired;
         case Double: return  ExplicitCastRequired;
         case List:   return  ExplicitCastRequired;
      }
   }

   SX_EXIT;
}


SxString SxVariant::getTypeName () const
{
   SX_TRACE ();
   switch (type) {
      case SxVariantType::Undefined:
         return "undefined";
      case SxVariantType::Any:
         return "any";
      case SxVariantType::Int:
         return "int";
      case SxVariantType::Double:
         return "float";
      case SxVariantType::Bool:
         return "bool";
      case SxVariantType::String:
         return "string";
      case SxVariantType::List:
         return "list";
      case SxVariantType::StringArray:
         return "StringArray";
      case SxVariantType::IntArray:
         return "IntArray";
      case SxVariantType::FloatArray:
         return "FloatArray";
      case SxVariantType::DoubleArray:
         return "DoubleArray";
      default: SX_EXIT;
   }
   return "unknown";
}

SxString SxVariant::getTypeArgs () const
{
   SX_TRACE ();
   SxString res;

   SxString typeStr;
   SxString rangeStr;
   SxString minStr;
   SxString maxStr;

   typeStr = getTypeName ();
   switch (type)  {
      case SxVariantType::Int:
         if (*(int64_t *)minPtr != SxVariantType::ScalarType<int>::minVal())
            minStr = SxString(*(int64_t *)minPtr);
         if (*(int64_t *)maxPtr != SxVariantType::ScalarType<int>::maxVal()) 
            maxStr = SxString(*(int64_t *)maxPtr);
         break;
      case SxVariantType::Double:
         if (*(double *)minPtr != SxVariantType::ScalarType<double>::minVal())
            minStr = SxString(*(double *)minPtr);
         if (*(double *)maxPtr != SxVariantType::ScalarType<double>::maxVal()) 
            maxStr = SxString(*(double *)maxPtr);
         break;
      case SxVariantType::String:
         if (minElem > 0)
            minStr = SxString(minElem);
         if (maxElem > 0)
            maxStr = SxString(maxElem);
         break;
      case SxVariantType::Bool:
         minStr = SxString ("false");
         maxStr = SxString ("true");
         break;
      case SxVariantType::List:
      case SxVariantType::Undefined:
      case SxVariantType::Any:
         break;
      default: SX_EXIT;
   }

   if (minStr != "" || maxStr != "")  {
      rangeStr = "(";
      if (minStr != "")
         rangeStr += "min " + minStr;
      if (maxStr != "")  {
         if (minStr != "")
            rangeStr += ", ";
         rangeStr += "max " + maxStr;
      }
      rangeStr += ")";
   }
   if (regex.getPtr() != NULL)
      rangeStr += "[/" + regexPattern + "/" + regexFlags + "]";

   return typeStr + rangeStr;
}

SxString SxVariant::getDescription (bool withTypeArgs) const
{
   SX_TRACE ();
   SxString res;

   if (tag != "")  {
      if (withTypeArgs)  res = tag + ": ";
      else               res = "<" + tag + ">";
   }

   SxString argsStr;
   if (withTypeArgs)  argsStr = getTypeArgs ();

   return res + argsStr;
}

SxString SxVariant::printToString () const
{
   SX_TRACE ();
   if (type != SxVariantType::List && !initialized)
      return "uninitialized";

   switch (type)  {
      case SxVariantType::Undefined:
         return "undefined";
      case SxVariantType::Any:
         return "any";
      case SxVariantType::Int:
         return SxString(*(int64_t *)data);
      case SxVariantType::Double:
         return SxString(*(double *)data);
      case SxVariantType::Bool:
         return ( (*(bool *)data)==true ?
                  SxString("true") : SxString("false") );
      case SxVariantType::String:
         return (*(SxString *)data);
      case SxVariantType::List:  {
         auto it = ((SxList<SxVariant> *)data)->begin ();
         SxList<SxString> res;
         while (it.isValid ())  {
            res.append (it->printToString ());
            ++it;
         }
         return "[" + SxString::join (res, ", ") + "]";
      }
      case SxVariantType::StringArray:  {
         auto it = ((SxArray<SxString> *)data)->begin ();
         SxList<SxString> res;
         while (it.isValid ())  {
            res.append (*it);
         }
         return "[" + SxString::join (res, ", ") + "]";
      }
      case SxVariantType::IntArray:  {
         auto it = ((SxArray<int64_t> *)data)->begin ();
         SxList<SxString> res;
         while (it.isValid ())  {
            res.append (SxString(*it));
         }
         return "[" + SxString::join (res, ", ") + "]";
      }
      case SxVariantType::FloatArray:  {
         auto it = ((SxArray<float> *)data)->begin ();
         SxList<SxString> res;
         while (it.isValid ())  {
            res.append (SxString(*it));
         }
         return "[" + SxString::join (res, ", ") + "]";
      }
      case SxVariantType::DoubleArray:  {
         auto it = ((SxArray<double> *)data)->begin ();
         SxList<SxString> res;
         while (it.isValid ())  {
            res.append (SxString(*it));
         }
         return "[" + SxString::join (res, ", ") + "]";
      }
      default: SX_EXIT;
   }

   return "unknown";
}

void SxVariant::undefTest (const SxVariant &a, const SxVariant &b) const
{
   SX_TRACE ();
   undefTest (a);
   undefTest (b);
}

void SxVariant::undefTest (const SxVariant &a) const
{
   SX_TRACE ();
   if (a.getType () == SxVariantType::Undefined)
      SX_THROW ("Undefined type '" + a.getTag () + "'");

   if (a.getType () == SxVariantType::Any)
      SX_THROW ("Undefined any type '" + a.getTag () + "'");

   if (!a.isInitialized ())
      SX_THROW ("Uninitialized value '" + a.getTag () + "'");
}

const SxArray<SxString> &SxVariant::toStringArray () const
{
   SX_TRACE ();
   if (type != SxVariantType::StringArray)  {
      SX_THROW ("toStringArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<SxString> *)data;
}

SxArray<SxString> &SxVariant::toStringArray ()
{
   SX_TRACE ();
   if (type != SxVariantType::StringArray)  {
      SX_THROW ("toStringArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<SxString> *)data;
}

const SxArray<int64_t> &SxVariant::toIntArray () const
{
   SX_TRACE ();
   if (type != SxVariantType::IntArray)  {
      SX_THROW ("toIntArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<int64_t> *)data;
}

SxArray<int64_t> &SxVariant::toIntArray ()
{
   SX_TRACE ();
   if (type != SxVariantType::IntArray)  {
      SX_THROW ("toIntArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<int64_t> *)data;
}

const SxArray<float> &SxVariant::toFloatArray () const
{
   SX_TRACE ();
   if (type != SxVariantType::FloatArray)  {
      SX_THROW ("toFloatArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<float> *)data;
}

SxArray<float> &SxVariant::toFloatArray ()
{
   SX_TRACE ();
   if (type != SxVariantType::FloatArray)  {
      SX_THROW ("toFloatArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<float> *)data;
}

const SxArray<double> &SxVariant::toDoubleArray () const
{
   SX_TRACE ();
   if (type != SxVariantType::DoubleArray)  {
      SX_THROW ("toDoubleArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<double> *)data;
}

SxArray<double> &SxVariant::toDoubleArray ()
{
   SX_TRACE ();
   if (type != SxVariantType::DoubleArray)  {
      SX_THROW ("toDoubleArray: Variant type is " + getTypeName ());
   }
   return *(SxArray<double> *)data;
}

int SxVariant::getTypeSize () const
{
   SX_TRACE ();
   int nBytes = 0;
   using namespace SxVariantType;
   switch (type)  {
      case Undefined:
         nBytes = 0;
         break;
      case String:
         nBytes = sizeof(char);
         break;
      case Int:
         nBytes = sizeof(int64_t);
         break;
      case Double:
         nBytes = sizeof(double);
         break;
      case StringArray:
         nBytes = sizeof(char);
         break;
      case IntArray:
         nBytes = sizeof(int64_t);
         break;
      case FloatArray:
         nBytes = sizeof(float);
         break;
      case DoubleArray:
         nBytes = sizeof(double);
         break;
      default:
         SX_EXIT;
   }
   return nBytes;
}

SxArray<char> SxVariant::pack () const
{
   SX_TRACE ();
   using namespace SxVariantType;

   SxArray<char> res;
   ssize_t pos = 0;
   int typeSize = getTypeSize ();

   ssize_t nElem = 0;
   uint8_t *src = NULL;
   if (type == String)  {
      SxString *p = (SxString *)data;
      src   = (uint8_t*)p->elements;
      nElem = p->getSize ();
   }  else if (  type == Int
            //|| type == Float
              || type == Double)
   {
      src   = (uint8_t*)data;
      nElem = 1;
   }  else if (type == StringArray)  {
      SxArray<SxString> *p = (SxArray<SxString> *)data;
      nElem = p->getSize ();
   }  else if (type == IntArray)  {
      SxArray<int64_t> *p = (SxArray<int64_t> *)data;
      src   = (uint8_t*)p->elements;
      nElem = p->getSize ();
   }  else if (type == FloatArray)  {
      SxArray<float> *p = (SxArray<float> *)data;
      src   = (uint8_t*)p->elements;
      nElem = p->getSize ();
   }  else if (type == DoubleArray)  {
      SxArray<double> *p = (SxArray<double> *)data;
      src   = (uint8_t*)p->elements;
      nElem = p->getSize ();
   }

   // --- total packed bytes
   ssize_t nLen = 0;
   ssize_t nBytes = 0;
   if (  type == String
      || type == IntArray
      || type == FloatArray
      || type == DoubleArray)
   {
      nLen = SxSerializer::getSize ((int64_t)nElem);// returns uint8_t
      nBytes = 1 + nLen + (typeSize * nElem);
   }  else if (type == StringArray)  {
      nLen = SxSerializer::getSize ((int64_t)nElem);
      nBytes = 1 + nLen;
      const SxArray<SxString> &a = toStringArray ();
      for (ssize_t i=0; i < nElem; ++i)  {
         ssize_t nElemStr = a(i).getSize ();
         ssize_t nLenStr = SxSerializer::getSize ((int64_t)nElemStr);
         nBytes += 1 + nLenStr + nElemStr;
      }
   }  else  {
      nBytes = typeSize;
   }
   nBytes += 2;// for type and typeSize
   //res.resize (nBytes + 1);
   res.resize (nBytes);
   uint8_t *buffer = (uint8_t*)res.elements;

   // --- type
   buffer[0] = static_cast<uint8_t>(type);
   buffer[1] = static_cast<uint8_t>(typeSize);
   pos += 2;

   // --- nElem
   if (  type == String
      || type == StringArray
      || type == IntArray
      || type == FloatArray
      || type == DoubleArray)
   {
      buffer[pos] = static_cast<uint8_t>(nLen);
      pos += 1;
      SxSerializer::pack ((int64_t)nElem, buffer + pos);
      pos += nLen;
   }

   // --- data
   if (typeSize == 1)  {
      if (type == StringArray)  {
         const SxArray<SxString> &a = toStringArray ();
         for (ssize_t i=0; i < nElem; ++i)  {
            ssize_t nElemStr = a(i).getSize ();
            ssize_t nLenStr = SxSerializer::getSize ((int64_t)nElemStr);
            buffer[pos] = static_cast<uint8_t>(nLenStr);
            pos += 1;
            SxSerializer::pack ((int64_t)nElemStr, buffer + pos);
            pos += nLen;
            memcpy (buffer + pos, a(i).elements, (size_t)nElemStr);
            pos += nElemStr;
         }
      }  else  {
         memcpy (buffer + pos, src, (size_t)nElem);
      }
   }  else if (typeSize == 4)  {
      SxEndian::hostToBig32 (src, buffer + pos, nElem * typeSize);
   }  else if (typeSize == 8)  {
      SxEndian::hostToBig64 (src, buffer + pos, nElem * typeSize);
   }  else if (typeSize != 0)  {
      SX_EXIT;
   }

   return res;
}

ssize_t SxVariant::unpack (const char *data_, ssize_t len)
{
   SX_TRACE ();
   using namespace SxVariantType;

   ssize_t pos = 0;
   if (len < 2)  {
      return -1;
   }
   SX_CHECK (data_);
   const uint8_t *buffer = (const uint8_t*)data_;

   // --- type
   int t = buffer[0];
   int typeSize = buffer[1];
   pos += 2;

   // --- nElem
   int64_t n = 1;
   if (  t == String
      || t == StringArray
      || t == IntArray
      || t == FloatArray
      || t == DoubleArray)
   {
      ssize_t ret = 0;
      ret = SxSerializer::unpack (data_ + pos, len - pos, &n);
      if (ret < 0) {
         return -1;
      }
      pos += ret;
   }  else if (t == Undefined)  {
      n = 0;
   }
   if (len < pos + (n * typeSize))  {
      return -1;
   }
   // --- data
   if (t == Undefined)  {
      setType (Undefined);

   }  else if (t == String)  {
      setType (String);
      *(SxString *)data = SxString (data_ + pos, n);

   }  else if (t == Int)  {
      uint64_t v = SxEndian::bigToHost64 (buffer + pos);
      int64_t i64 = 0;
      memcpy (&i64, &v, sizeof v);
      setType (Int);
      *(int64_t *)data = i64;

   /*} else if (t == Float)  {
      uint32_t v = SxEndian::bigToHost32 (buffer + pos);
      float f = 0;
      memcpy (&f, &v, sizeof v);
      setType (Float);
      *(float *)data = f;*/

   }  else if (t == Double)  {
      uint64_t v = SxEndian::bigToHost64 (buffer + pos);
      double d = 0;
      memcpy (&d, &v, sizeof v);
      setType (Double);
      *(double *)data = d;

   }  else if (t == StringArray)  {
      setType (StringArray);
      SxArray<SxString> &a = toStringArray ();
      a.resize (n);
      for (ssize_t i=0; i < n; ++i)  {
         ssize_t r = SxSerializer::unpack (data_ + pos, len - pos, &a(i));
         if (r < 0)  {
            return -1;
         }
         pos += r;
      }

   }  else if (t == IntArray)  {
      ssize_t nBytes = n * typeSize;
      SxArray<uint8_t> mem(nBytes);
#     ifdef HAS_BIG_ENDIAN
         memcpy (mem.elements, buffer + pos, (size_t)nBytes);
#     else
         SxEndian::swap64 (buffer + pos, mem.elements, nBytes);
#     endif /* HAS_BIG_ENDIAN */
      setType (IntArray);
      SxArray<int64_t> *p = (SxArray<int64_t> *)data;
      *p = SxArray<int64_t>(reinterpret_cast<int64_t*>(mem.elements), n);

   }  else if (t == FloatArray)  {
      ssize_t nBytes = n * typeSize;
      SxArray<uint8_t> mem(nBytes);
#     ifdef HAS_BIG_ENDIAN
         memcpy (mem.elements, buffer + pos, (size_t)nBytes);
#     else
         SxEndian::swap32 (buffer + pos, mem.elements, nBytes);
#     endif /* HAS_BIG_ENDIAN */
      setType (FloatArray);
      SxArray<float> *p = (SxArray<float> *)data;
      *p = SxArray<float>(reinterpret_cast<float*>(mem.elements), n);

   }  else if (t == DoubleArray)  {
      ssize_t nBytes = n * typeSize;
      SxArray<uint8_t> mem(nBytes);
#     ifdef HAS_BIG_ENDIAN
         memcpy (mem.elements, buffer + pos, (size_t)nBytes);
#     else
         SxEndian::swap64 (buffer + pos, mem.elements, nBytes);
#     endif /* HAS_BIG_ENDIAN */
      setType (DoubleArray);
      SxArray<double> *p = (SxArray<double> *)data;
      *p = SxArray<double>(reinterpret_cast<double*>(mem.elements), n);
   }
   if (t != StringArray)  {
      pos += n * typeSize;
   }

   return pos;
}

// --------------------------------------------------------------------------

SxString SxVariant::getTypeStr (int type)
{
   switch (type)  {
      case SxVariantType::Undefined:
         return "undefined";
      case SxVariantType::Any:
         return "any";
      case SxVariantType::Int:
         return "int";
      case SxVariantType::Double:
         return "float";
      case SxVariantType::Bool:
         return "bool";
      case SxVariantType::String:
         return "string";
      case SxVariantType::List:
         return "list";
      case SxVariantType::StringArray:
         return "StringArray";
      case SxVariantType::IntArray:
         return "IntArray";
      case SxVariantType::FloatArray:
         return "FloatArray";
      case SxVariantType::DoubleArray:
         return "DoubleArray";
      default: SX_EXIT;
   }
   return "unknown";
}

int SxVariant::getTypeId (const SxString &str)
{
   if (str == "list")  {
      return SxVariantType::List;
   }  else if (str == "int")  {
      return SxVariantType::Int;
   }  else if (str == "float")  {
      return SxVariantType::Double;
   }  else if (str == "string")  {
      return SxVariantType::String;
   }  else if (str == "bool")  {
      return SxVariantType::Bool;
   }  else if (str == "StringArray")  {
      return SxVariantType::StringArray;
   }  else if (str == "IntArray")  {
      return SxVariantType::IntArray;
   }  else if (str == "FloatArray")  {
      return SxVariantType::FloatArray;
   }  else if (str == "DoubleArray")  {
      return SxVariantType::DoubleArray;
   }  else  {
      return SxVariantType::Undefined;
   }
}

#define SX_VARIANT_SET_FUNC_IMPLEMENTATION(T)                                  \
   void SxVariant::set (const T in)                                            \
   {                                                                           \
      SX_CHECK (dataType == SxVariantType::ScalarType<T>::Type);               \
                                                                               \
      typedef typename SxVariantType::ScalarType<T>::AtomType A;               \
      if (in < *(A *)minPtr || in > *(A *)maxPtr)                              \
         SX_THROW("Value " + SxString(in) + " is out of range");               \
                                                                               \
      *(A *)data = in;                                                         \
      initialized = true;                                                      \
   } \
   struct SX_UNIQUE_ID(whatever) { }

SX_VARIANT_SET_FUNC_IMPLEMENTATION(int32_t);
SX_VARIANT_SET_FUNC_IMPLEMENTATION(int64_t);
// disable warning for unsave use of [<,>] against bool type
#ifdef MSVC
#  pragma warning (suppress:4804)
   SX_VARIANT_SET_FUNC_IMPLEMENTATION(bool);
#else
   SX_VARIANT_SET_FUNC_IMPLEMENTATION(bool);
#endif

SX_VARIANT_SET_FUNC_IMPLEMENTATION(float);
SX_VARIANT_SET_FUNC_IMPLEMENTATION(double);

// --- List
void SxVariant::append (const SxVariant &in)
{
   SX_CHECK (dataType == SxVariantType::List);

   if ((*(int *)minPtr && in.getType () < *(int *)minPtr)
       || (*(int *)maxPtr && in.getType () > *(int *)maxPtr))
      SX_THROW("The value to append does not match list type");

   ((SxList<SxVariant> *)data)->append (in);
}

