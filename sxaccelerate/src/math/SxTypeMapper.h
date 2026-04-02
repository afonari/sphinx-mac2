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

#ifndef _SX_TYPE_MAPPER_H_
#define _SX_TYPE_MAPPER_H_
#include <SxComplex.h>

template<class T> class SxTypeMapper;
// --- The next 60 lines were generated from snippets/SxVector.h snippet SxTypeMapper
// --- template specialization for float
template<>
class SxTypeMapper<float> {
   public:
      /// The real (non-complex) scalar type
      typedef float TReal;
      /// The complex scalar type
      typedef SxComplex8 TCmplx;
};

// --- template specialization for double
template<>
class SxTypeMapper<double> {
   public:
      /// The real (non-complex) scalar type
      typedef double TReal;
      /// The complex scalar type
      typedef SxComplex16 TCmplx;
};

// --- template specialization for SxComplex8
template<>
class SxTypeMapper<SxComplex8> {
   public:
      /// The real (non-complex) scalar type
      typedef float TReal;
      /// The complex scalar type
      typedef SxComplex8 TCmplx;
};

// --- template specialization for SxComplex16
template<>
class SxTypeMapper<SxComplex16> {
   public:
      /// The real (non-complex) scalar type
      typedef double TReal;
      /// The complex scalar type
      typedef SxComplex16 TCmplx;
};

// --- template specialization for int
template<>
class SxTypeMapper<int> {
   public:
      /// The real (non-complex) scalar type
      typedef int TReal;
      /// The complex scalar type
      typedef SxComplex<int> TCmplx;
};

// --- template specialization for long
template<>
class SxTypeMapper<long> {
   public:
      /// The real (non-complex) scalar type
      typedef long TReal;
      /// The complex scalar type
      typedef SxComplex<long> TCmplx;
};
// --- SxTypeMapper

#endif
