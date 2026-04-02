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


SxVariant::SxVariant ()
   : type(SxVariantType::Undefined),
     dataType(SxVariantType::Undefined),
     initialized(false),
     data(NULL),
     minPtr(NULL),
     maxPtr(NULL),
     minElem(0),
     maxElem(0)
{
   // empty
}

template<class T>
SxVariant::SxVariant (const T &in, const SxString &tag_)
   : type(SxVariantType::Undefined),
     dataType(SxVariantType::Undefined),
     initialized(false),
     data(NULL),
     minPtr(NULL),
     maxPtr(NULL),
     minElem(0),
     maxElem(0),
     tag(tag_)
{
   setTypeT (in);
   set (in);
}

template<class A, class B, class T>
SxVariant::SxVariant (const T &in,
                      const A &min,
                      const B &max,
                      const SxString &tag_)
   : type(SxVariantType::Undefined),
     dataType(SxVariantType::Undefined),
     initialized(false),
     data(NULL),
     minPtr(NULL),
     maxPtr(NULL),
     minElem(0),
     maxElem(0),
     tag(tag_)
{
   setTypeT (in);
   setLimits (min, max);
   set (in);
}

SxVariant::SxVariant (const SxVariant &in)
   : type(SxVariantType::Undefined),
     dataType(SxVariantType::Undefined),
     initialized(false),
     data(NULL),
     minPtr(NULL),
     maxPtr(NULL),
     minElem(0),
     maxElem(0)
{
   *this = in;
}

// --- arrays
template<class T>
SxVariant::SxVariant (const SxArray<T> &in)
   : type(SxVariantType::Undefined),
     data(NULL)
{
   set (in);
}

#ifdef MSVC
   // Constructor in Constructor is not recognized as initialization -> disable warning
#  pragma warning( push )
#  pragma warning( disable : 26495)
#endif

SxVariant::SxVariant (int32_t in)
   : SxVariant (in, "")
{
   // empty
}

SxVariant::SxVariant (int64_t in)
   : SxVariant (in, "")
{
   // empty
}

SxVariant::SxVariant (bool in)
   : SxVariant (in, "")
{
   // empty
}

SxVariant::SxVariant (float in)
   : SxVariant (in, "")
{
   // empty
}

SxVariant::SxVariant (double in)
   : SxVariant (in, "")
{
   // empty
}

SxVariant::SxVariant (const char *in)
   : SxVariant (in, "")
{
   // empty
}

SxVariant::SxVariant (const SxString &in)
   : SxVariant (in, "")
{
   // empty
}

SxVariant::SxVariant (const SxList<SxVariant> &in)
   : SxVariant (in, "")
{
   // empty
}

#ifdef MSVC
#  pragma warning( pop )
#endif

SxVariant::~SxVariant ()
{
   reset ();
}

bool SxVariant::isInitialized () const
{
   return initialized;
}

int SxVariant::getType () const
{
   return type;
}

int SxVariant::getDataType () const
{
   return dataType;
}

template<class T>
void SxVariant::setTypeT (const T)
{
   setType (SxVariantType::ScalarType<T>::Type);
}

template<>
inline void SxVariant::setTypeT (const SxVariantType::DataType in)
{
   setType (in);
}

template<>
inline void SxVariant::setTypeT (const char *)
{
   setType (SxVariantType::String);
}

void SxVariant::setTypeT (const SxString &)
{
   setType (SxVariantType::String);
}

// --- Limits
template<class A, class B>
SxVariant &SxVariant::setLimits (const A &a, const B &b)
{
   SX_CHECK (dataType != SxVariantType::Undefined);

   if (dataType == SxVariantType::String)  {
      minElem = (ssize_t)a;
      maxElem = (ssize_t)b;
      return *this;
   }

   switch (dataType) {
      case SxVariantType::Int:
         *(int64_t *)minPtr = (int64_t)a;
         *(int64_t *)maxPtr = (int64_t)b;
         break;
      case SxVariantType::Double:
         *(double *)minPtr = (double)a;
         *(double *)maxPtr = (double)b;
         break;
      case SxVariantType::String:
         *(char *)minPtr = (char)a;
         *(char *)maxPtr = (char)b;
         break;
      case SxVariantType::List:
         *(int *)minPtr = (int)a;
         *(int *)maxPtr = (int)b;
         break;
      default: SX_EXIT;
   }
   return *this;
}

template<class A, class B>
SxVariant &SxVariant::getLimits (A *a, B *b)
{
   SX_CHECK (dataType != SxVariantType::Undefined);

   switch (dataType) {
      case SxVariantType::Int:
         *a = (A)*(int64_t *)minPtr;
         *b = (B)*(int64_t *)maxPtr;
         break;
      case SxVariantType::Double:
         *a = (A)*(double *)minPtr;
         *b = (B)*(double *)maxPtr;
         break;
      case SxVariantType::String:
         *a = (A)minElem;
         *b = (B)maxElem;
         break;
      case SxVariantType::List:
         *a = (A)*(int *)minPtr;
         *b = (B)*(int *)maxPtr;
         break;
      default: SX_EXIT;
   }
   return *this;
}

template<class A>
SxVariant &SxVariant::setMin (const A &min)
{
   SX_CHECK (dataType != SxVariantType::Undefined);

   if (dataType == SxVariantType::String)  {
      minElem = (ssize_t)min;
      return *this;
   }

   switch (dataType) {
      case SxVariantType::Int:
         *(int64_t *)minPtr = (int64_t)min;
         break;
      case SxVariantType::Double:
         *(double *)minPtr = (double)min;
         break;
      case SxVariantType::String:
         *(char *)minPtr = (char)min;
         break;
      case SxVariantType::List:
         *(int *)minPtr = (int)min;
         break;
      default: SX_EXIT;
   }
   return *this;
}

template<class A>
SxVariant &SxVariant::setMax (const A &max)
{
   SX_CHECK (dataType != SxVariantType::Undefined);

   if (dataType == SxVariantType::String)  {
      maxElem = (ssize_t)max;
      return *this;
   }

   switch (dataType) {
      case SxVariantType::Int:
         *(int64_t *)maxPtr = (int64_t)max;
         break;
      case SxVariantType::Double:
         *(double *)maxPtr = (double)max;
         break;
      case SxVariantType::String:
         *(char *)maxPtr = (char)max;
         break;
      case SxVariantType::List:
         *(int *)maxPtr = (int)max;
         break;
      default: SX_EXIT;
   }
   return *this;
}

SxVariant &SxVariant::setRegex (const SxString &regex_,
                                const SxString &flags_)
{
   regex = SxPtr<SxRegex>::create (regex_, flags_);
   regexPattern = regex_;
   regexFlags   = flags_;
   return *this;
}

SxString SxVariant::getRegex () const
{
   return "/" + regexPattern + "/" + regexFlags;
}

void SxVariant::set (const SxVariantType::DataType in)
{
   setTypeT (in);
}

template<>
inline void SxVariant::set (const SxArray<SxString> &in)
{
   setType (SxVariantType::StringArray);
   *(SxArray<SxString> *)data = in;
   initialized = true;
}

template<class T>
void SxVariant::set (const SxArray<T> &in)
{
   setType (SxVariantType::ArrayType<T>::Type);

   T *srcPtr = in.elements;
   ssize_t nElem = in.getSize ();

   if (type == SxVariantType::IntArray) {
      SxArray<int64_t> *p = (SxArray<int64_t> *)data;
      p->resize (nElem);
      int64_t *dstPtr = p->elements;
      for (ssize_t i=0; i < nElem; ++i) {
         *dstPtr++ = static_cast<int64_t>(*srcPtr++);
      }
      initialized = true;
   } else if (type == SxVariantType::FloatArray) {
      SxArray<float> *p = (SxArray<float> *)data;
      p->resize (nElem);
      float *dstPtr = p->elements;
      for (ssize_t i=0; i < nElem; ++i) {
         *dstPtr++ = static_cast<float>(*srcPtr++);
      }
      initialized = true;
   } else if (type == SxVariantType::DoubleArray) {
      SxArray<double> *p = (SxArray<double> *)data;
      p->resize (nElem);
      double *dstPtr = p->elements;
      for (ssize_t i=0; i < nElem; ++i) {
         *dstPtr++ = static_cast<double>(*srcPtr++);
      }
      initialized = true;
   }
}

void SxVariant::set (const char *in)
{
   SX_CHECK (dataType == SxVariantType::String);

   if (hasStringLimits ())
      validate (in);

   *(SxString *)data = in;
   initialized = true;
}

void SxVariant::set (const SxString &in)
{
   SX_CHECK (dataType == SxVariantType::String);

   if (hasStringLimits ())
      validate (in);

   *(SxString *)data = in;
   initialized = true;
}

void SxVariant::set (const SxVariant &in, bool allowTypecast)
{
   if (in.dataType != dataType && allowTypecast)  {
      switch (dataType) {
         case SxVariantType::Int:
            set (in.toInt ());
            break;
         case SxVariantType::Double:
            set (in.toDouble ());
            break;
         case SxVariantType::String:
            set (in.toString ());
            break;
         default: SX_EXIT;
      }
   } else {
      switch (in.dataType) {
         case SxVariantType::Undefined:
            // TODO:
            break;
         case SxVariantType::Int:
            set (in.getInt ());
            break;
         case SxVariantType::Double:
            set (in.getDouble ());
            break;
         case SxVariantType::String:
            set (in.getString ());
            break;
         case SxVariantType::List:
            set (*(SxList<SxVariant> *)in.data);
            break;
         case SxVariantType::StringArray:
            *(SxArray<SxString> *)data = in.toStringArray ();
            initialized = true;
            break;
         case SxVariantType::IntArray:
            *(SxArray<int64_t> *)data = in.toIntArray ();
            initialized = true;
            break;
         case SxVariantType::FloatArray:
            *(SxArray<float> *)data = in.toFloatArray ();
            initialized = true;
            break;
         case SxVariantType::DoubleArray:
            *(SxArray<double> *)data = in.toDoubleArray ();
            initialized = true;
            break;
         default: SX_EXIT;
      }
   }
   //tag = in.tag;
}

void SxVariant::set (const SxList<SxVariant> &in)
{
   SX_CHECK (dataType == SxVariantType::List);

   SxList<SxVariant>::ConstIterator it = in.begin ();
   while (it.isValid ())  {
      append (*it);
      ++it;
   }

   initialized = true;
}

template<class T>
void SxVariant::operator= (const T in)
{
   setType (SxVariantType::ScalarType<T>::Type);
   set (in);
}

void SxVariant::operator= (const SxString &in)
{
   setType (SxVariantType::String);
   set (in);
}

template<class T>
void SxVariant::operator= (const SxArray<T> &in)
{
   setType (SxVariantType::ArrayType<T>::Type);
   set (in);
}

int SxVariant::getRank () const
{
   if (dataType == SxVariantType::List)  { // [....]
      if (getListSize () == 0)  {          // []
         return 1;  // empty vector
      }  else if (getListSize () == 1)  {  // [[...]]
         return 1 + begin()->getRank();
      }  else  {
         SxList<SxVariant>::ConstIterator it = begin ();
         int rk, maxRank = 0;
         while (it.isValid ())  {
            rk = (*it).getRank();
            if (rk > maxRank)  maxRank = rk;
            ++it;
         }
         return 1 + maxRank;
      }
      SX_EXIT;
      return -1;
   }  else  {
      return 0;                               // scalar 123.456;
   }
   return 0;                                  // scalar
}

ssize_t SxVariant::getListSize () const
{
   return (  dataType == SxVariantType::List)
           ? ((SxList<SxVariant> *)data)->getSize ()
           : 0;
}

SxList<SxVariant>::Iterator SxVariant::begin ()
{
   return (  dataType == SxVariantType::List)
           ? ((SxList<SxVariant> *)data)->begin ()
           : SxList<SxVariant>::Iterator ();
}

SxList<SxVariant>::ConstIterator SxVariant::begin () const
{
   return (  dataType == SxVariantType::List)
           ? ((SxList<SxVariant> *)data)->begin ()
           : SxList<SxVariant>::Iterator ();
}

SxList<SxVariant>::Iterator SxVariant::end ()
{
   return (  dataType == SxVariantType::List)
           ? ((SxList<SxVariant> *)data)->end ()
           : SxList<SxVariant>::Iterator ();
}

SxList<SxVariant>::ConstIterator SxVariant::end () const
{
   return (  dataType == SxVariantType::List)
           ? ((SxList<SxVariant> *)data)->end ()
           : SxList<SxVariant>::Iterator ();
}

int64_t SxVariant::getInt () const
{
   SX_CHECK (dataType == SxVariantType::Int, getTypeName ());
   return *(int64_t *)data;
}

double SxVariant::getDouble () const
{
   SX_CHECK (dataType == SxVariantType::Double, getTypeName ());
   return *(double *)data;
}

bool SxVariant::getBool () const
{
   SX_CHECK (dataType == SxVariantType::Bool, getTypeName ());
   return *(bool *)data;
}

const SxString &SxVariant::getString () const
{
   SX_CHECK (dataType == SxVariantType::String, getTypeName ());
   return *(SxString *)data;
}

int64_t SxVariant::toInt () const
{
   if (dataType == SxVariantType::Int)
      return *(int64_t *)data;
   else if (dataType == SxVariantType::Double)
      return (int64_t)lround(*(double *)data);
   else if (dataType == SxVariantType::String)
      return ((SxString *)data)->toInt ();
   else
      SX_CHECK (false, "typecast not possible");

   return  0;
}

double SxVariant::toDouble () const
{
   if (dataType == SxVariantType::Int)
      return (double)(*(int64_t *)data);
   else if (dataType == SxVariantType::Double)
      return *(double *)data;
   else if (dataType == SxVariantType::String)
      return ((SxString *)data)->toDouble ();
   else
      SX_CHECK (false, "typecast not possible");

   return  0.0;
}

bool SxVariant::toBool () const
{
   if (dataType == SxVariantType::Int)
      return (bool)(*(int64_t *)data);
   else if (dataType == SxVariantType::Double)
      return (bool)(*(double *)data);
   else if (dataType == SxVariantType::Bool)
      return *(bool *)data;
   else
      SX_CHECK (false, "typecast not possible");

   return false;
}

SxString SxVariant::toString () const
{
   if (dataType == SxVariantType::Int)
      return SxString(*(int64_t *)data);
   else if (dataType == SxVariantType::Double)
      return SxString(*(double *)data);
   else if (dataType == SxVariantType::String)
      return *(SxString *)data;
   else
      SX_CHECK (false, "typecast not possible");

   return SxString();
}

SxString SxVariant::toByteArray () const
{
   switch (dataType)  {
      case SxVariantType::Int:
         return SxString(*(int64_t *)data);
      case SxVariantType::Double:
         return SxString(*(double *)data);
      case SxVariantType::String:
         return *(SxString *)data;
      default: SX_EXIT;
   }

   return SxString();
}

int64_t SxVariant::minInt () const
{
   return *(int64_t *)minPtr;
}

int64_t SxVariant::maxInt () const
{
   return *(int64_t *)maxPtr;
}

double SxVariant::minDouble () const
{
   return *(double *)minPtr;
}

double SxVariant::maxDouble () const
{
   return *(double *)maxPtr;
}

const SxString &SxVariant::getTag () const
{
   return tag;
}

void SxVariant::setTag (const SxString &tag_)
{
   tag = tag_;
}
#ifdef WIN32
std::wostream& operator<< (std::wostream &s, const SxVariant &in)
{
   s << in.getDescription (true) << " " << in.printToString ();
   return s;
}
#endif
std::ostream& operator<< (std::ostream &s, const SxVariant &in)
{
   s << in.getDescription (true) << " " << in.printToString ();
   return s;
}

bool SxVariant::operator== (const SxVariant &b) const
{
   const SxVariant &a = *this;
   int typeA = a.getDataType(), typeB = b.getDataType();

   if (   typeA == SxVariantType::Int
       && typeB == SxVariantType::Int )  {
      return a.getInt() == b.getInt();
   } else if (   typeA == SxVariantType::Bool
              && typeB == SxVariantType::Bool )  {
      return a.getBool() == b.getBool();
   } else if (   typeA == SxVariantType::String
              && typeB == SxVariantType::String )  {
      return a.getString() == b.getString();
   } else if (   typeA == SxVariantType::Undefined
              && typeB == SxVariantType::Undefined )  {
      return true;
   } else if (   (   typeA == SxVariantType::Undefined
                  && typeB != SxVariantType::Undefined)
              || (   typeA != SxVariantType::Undefined
                  && typeB == SxVariantType::Undefined) )  {
      return false;
   }

   SX_EXIT; // no typecast available
   return false;
}

bool SxVariant::operator!= (const SxVariant &b) const
{
   return !operator== (b);
}

bool SxVariant::operator< (const SxVariant &b) const
{
   const SxVariant &a = *this;
   int typeA = a.getType(), typeB = b.getType();

   if (  (typeA == SxVariantType::Int || typeA == SxVariantType::Double)
      && (typeB == SxVariantType::Int || typeB == SxVariantType::Double))
   {
      if (typeA == SxVariantType::Double || typeB == SxVariantType::Double)
         return a.toDouble() < b.toDouble();
      else
         return a.toInt() < b.toInt();
   }

   SX_EXIT; // no typecast available
   return false;
}
