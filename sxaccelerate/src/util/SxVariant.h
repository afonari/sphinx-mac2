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

#ifndef _SX_VARIANT_H_
#define _SX_VARIANT_H_

#include <SxUtil.h>
#include <SxPtr.h>
#include <SxRegex.h>

namespace SxVariantType {

   enum DataType { Undefined, Any, Int, Double, Bool, String, List,
                   StringArray, IntArray, FloatArray, DoubleArray };

   template<class T, size_t N=0> class ScalarType {
      public: enum { Type = Undefined };
   };

   template<> class ScalarType<int32_t> {
      public:
         enum { Type = Int };
         typedef int64_t AtomType;
         inline static int64_t minVal () { return INT64_MIN ;}
         inline static int64_t maxVal () { return INT64_MAX; }
   };
   template<> class ScalarType<int64_t> {
      public:
         enum { Type = Int };
         typedef int64_t AtomType;
         inline static int64_t minVal () { return INT64_MIN ;}
         inline static int64_t maxVal () { return INT64_MAX; }
   };
   template<> class ScalarType<bool> {
      public:
         enum { Type = Bool };
         typedef uint8_t AtomType;
         inline static uint8_t minVal () { return 0; }
         inline static uint8_t maxVal () { return 1; }
   };

   template<> class ScalarType<float> {
      public:
         enum { Type = Double };
         typedef double AtomType;
         inline static float minVal () { return -HUGE_VAL; }
         inline static float maxVal () { return  HUGE_VAL; }
   };
   template<> class ScalarType<double> {
      public:
         enum { Type = Double };
         typedef double AtomType;
         inline static double minVal () { return -HUGE_VAL; }
         inline static double maxVal () { return  HUGE_VAL; }
   };

   template<> class ScalarType<const char *> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };
   template<> class ScalarType<char[]> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };
   template<size_t N> class ScalarType<char[N]> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };
   template<> class ScalarType<SxString> {
      public:
         enum { Type = String };
         typedef char AtomType;
         inline static char minVal () { return CHAR_MIN; }
         inline static char maxVal () { return CHAR_MAX; }
   };

   template<class T> class ArrayType {
      public: enum { Type = Undefined };
   };
   template<> class ArrayType<SxString> : public ScalarType<SxString> {
      public:
         enum { Type = StringArray };
         typedef SxArray<SxString> Array;
   };
   template<> class ArrayType<int32_t> : public ScalarType<int32_t> {
      public:
         enum { Type = IntArray };
         typedef SxArray<int64_t> Array;
   };
   template<> class ArrayType<int64_t> : public ScalarType<int64_t> {
      public:
         enum { Type = IntArray };
         typedef SxArray<int64_t> Array;
   };
   template<> class ArrayType<float> : public ScalarType<float> {
      public:
         enum { Type = FloatArray };
         typedef SxArray<float> Array;
   };
   template<> class ArrayType<double> : public ScalarType<double> {
      public:
         enum { Type = DoubleArray };
         typedef SxArray<double> Array;
   };
   template<> class ArrayType<char> : public ScalarType<char []> {
      public:
         enum { Type = String };
   };
}

/** \brief Variant of simplest data types

    This class provides a variant of simple data types. Using this variant
    lists or arrays with elements of varying types can be created at runtime.

    \author Sixten Boeck */
class SX_EXPORT_UTIL SxVariant
{
   public:

      enum CastType { ExactCast, ImplicitCast, ExplicitCastRequired,
                      IncompatibleTypes };

      inline SxVariant ();

      inline SxVariant (const SxVariant &in);
      void operator= (const SxVariant &in);

      template<class T>
      inline explicit SxVariant (const T &in, const SxString &tag="");

      template<class A, class B, class T>
      inline  SxVariant (const T &in,
                         const A &min,
                         const B &max,
                         const SxString &tag="");

      // --- arrays
      template<class T> inline SxVariant (const SxArray<T> &in);
      template<class T> inline void operator= (const SxArray<T> &in);
      template<class T> inline void set (const SxArray<T> &in);

      inline SxVariant (int32_t in);
      inline SxVariant (int64_t in);
      inline SxVariant (bool in);
      inline SxVariant (float in);
      inline SxVariant (double in);
      inline SxVariant (const char *in);
      inline SxVariant (const SxString &in);
      inline SxVariant (const SxList<SxVariant> &in);

      inline ~SxVariant ();

      inline bool  operator== (const SxVariant &) const;
      inline bool  operator!= (const SxVariant &) const;
      inline bool  operator< (const SxVariant &) const;

      // --- Set value (with type/range/regex check)
      template<class T> inline void operator= (const T in);
                        inline void operator= (const SxString &in);

      // new set functions
      void set (int32_t);
      void set (int64_t);
      void set (bool);
      void set (float);
      void set (double);
      inline void set (const char *in);
      inline void set (const SxString &in);
      inline void set (const SxVariant &, bool allowTypecast=false);
      inline void set (const SxList<SxVariant> &in);
      inline void set (const SxVariantType::DataType in);

      // --- List
      void append (const SxVariant &in);

      inline int     getRank () const;
      inline ssize_t getListSize () const;

      // --- List iterators
      inline SxList<SxVariant>::Iterator begin ();
      inline SxList<SxVariant>::ConstIterator begin () const;
      inline SxList<SxVariant>::Iterator end ();
      inline SxList<SxVariant>::ConstIterator end () const;

      // --- Limits
      template<class A, class B>
      inline SxVariant &setLimits (const A &min, const B &max);

      template<class A, class B>
      inline SxVariant &getLimits (A *min, B *max);

      template<class A> inline SxVariant &setMin (const A &min);
      template<class A> inline SxVariant &setMax (const A &max);

      void resetLimits ();

      inline SxVariant &setRegex (const SxString &, const SxString & = "");
      inline SxString   getRegex () const;

      // --- Match
      bool match (const SxVariant &in, bool withLimits=false) const;
      bool matchLimits (const SxVariant &in) const;
      CastType castDiff (const SxVariant &in) const;

      // --- Values
      inline int64_t getInt () const;
      inline double getDouble () const;
      inline bool   getBool () const;
      inline const SxString &getString () const;

      // --- Values with cast
      inline int64_t toInt () const;
      inline double toDouble () const;
      inline bool   toBool ()  const;
      inline SxString toString () const;
      inline SxString toByteArray () const;

      // --- arrays
      const SxArray<SxString> &toStringArray () const;
      SxArray<SxString> &toStringArray ();

      const SxArray<int64_t> &toIntArray () const;
      SxArray<int64_t> &toIntArray ();

      const SxArray<float> &toFloatArray () const;
      SxArray<float> &toFloatArray ();

      const SxArray<double> &toDoubleArray () const;
      SxArray<double> &toDoubleArray ();

      int getTypeSize () const;

      // --- explicit type initialization
      void setType (int type);

      inline bool isInitialized () const;
      inline int  getType () const;
      inline int  getDataType () const;
      inline const SxString &getTag () const;
      inline void setTag (const SxString &);

      SxString printToString () const;
      SxString getTypeArgs () const;
      SxString getTypeName () const;
      SxString getDescription (bool withTypeArgs=false) const;

      SxArray<char> pack () const;
      ssize_t unpack (const char *data, ssize_t len);


      static SxString getTypeStr (int type);
      static int      getTypeId  (const SxString &str);

   protected:

      int type;
      int dataType;
      bool initialized;
      void *data;

      void *minPtr, *maxPtr;
      ssize_t minElem;
      ssize_t maxElem;
      SxPtr<SxRegex> regex;
      SxString regexPattern, regexFlags;

      SxString tag;

      template<class T> inline void setTypeT (const T in);
                        inline void setTypeT (const SxString &in);


      bool hasStringLimits () const;
      bool isValidString (const SxString &) const;
      void validate (const SxString &) const;

      void undefTest (const SxVariant &a, const SxVariant &b) const;
      void undefTest (const SxVariant &a) const;

      void reset ();

      inline int64_t minInt () const;
      inline int64_t maxInt () const;
      inline double minDouble () const;
      inline double maxDouble () const;
};

#ifdef WIN32
   SX_EXPORT_UTIL inline std::wostream &operator<< (std::wostream &,
                                                    const SxVariant &);
#endif
SX_EXPORT_UTIL inline std::ostream &operator<< (std::ostream &,
                                                const SxVariant &);

#include <SxVariant.hpp>

#endif /* _SX_VARIANT_H_ */
