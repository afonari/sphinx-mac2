#ifndef SXVEC_INSTANCE
#define SXVEC_INSTANCE extern template
#endif

SXVEC_INSTANCE typename SxVecResult<double,Compact>::VecType SxVecCompute<Simd>::abs (const SxVecRef<SxComplex<double>, Compact> &vec);
SXVEC_INSTANCE typename SxVecResult<double,Compact>::VecType SxVecCompute<Simd>::absSqr (const SxVecRef<SxComplex<double>, Compact> &vec);

// Generate a list of instantiation candidates from SxVecCompute.hpp
// sed -ne 's/.* \(\S*::\S*\) .*/\1/p' SxVecCompute_instantiate.h | sort | uniq > instantiate_list
// sed -ne'/^template/{s/template *<[^>]*>/SXVEC_INSTANCE/;:loop;N;s/{.*/;/;tend;bloop;:end;s/\n */ /g;p;}' SxVecCompute.hpp | grep -v inline | grep -vf instantiate_list
//
/** --- Explicit instantiations
    Generate a complete list via (bash)
 { for id in `./snippets/codegen.pl snippets/SxVecCompute_instantiate 2>&1 | grep -v :` ; do ./snippets/codegen.pl snippets/SxVecCompute_instantiate --id $id ; done ; }
    ---
 */

// --- The next 73 lines were generated from snippets/SxVecCompute_instantiate snippet UseIt
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, SubMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, Strided>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, Strided>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, Strided>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, Strided>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, Strided>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, Strided>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, Strided>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, Strided>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, Strided>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, Strided>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, Strided>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, Strided>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, SubMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, Strided> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const SxVecRef<double, GeneralMatrix>  &x, const SxVecRef<double, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, Strided>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, SubMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, Strided> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxVecRef<SxComplex16, GeneralMatrix> &y);
// --- UseIt
// --- The next 49 lines were generated from snippets/SxVecCompute_instantiate snippet UseIt_scalar
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::Divide<true>::divide (const SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::Divide<true>::divideInPlace (SxVecRef<double, SubMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::Divide<true>::divide (const SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::Divide<true>::divideInPlace (SxVecRef<double, Strided>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::add (const SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::multiply (const SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::Divide<true>::divide (const SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::Divide<true>::divideInPlace (SxVecRef<double, GeneralMatrix>  &x, const double &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::Divide<true>::divide (const SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::Divide<true>::divideInPlace (SxVecRef<SxComplex16, SubMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::Divide<true>::divide (const SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::Divide<true>::divideInPlace (SxVecRef<SxComplex16, Strided>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::add (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::multiply (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::Divide<true>::divide (const SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::Divide<true>::divideInPlace (SxVecRef<SxComplex16, GeneralMatrix>  &x, const SxComplex16 &y) ;
// --- UseIt_scalar
// --- The next 19 lines were generated from snippets/SxVecCompute_instantiate snippet Simd_scalar_subdiv
SXVEC_INSTANCE SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::subtract (const double &y, const SxVecRef<double, Compact> &x) ;
SXVEC_INSTANCE SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::divide (const double &y, const SxVecRef<double, Compact> &x) ;
SXVEC_INSTANCE SxVecResult<double, PackedSymMatrix>::VecType SxVecCompute<Simd>::subtract (const double &y, const SxVecRef<double, PackedSymMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<double, PackedSymMatrix>::VecType SxVecCompute<Simd>::divide (const double &y, const SxVecRef<double, PackedSymMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<double, BandedMatrix>::VecType SxVecCompute<Simd>::subtract (const double &y, const SxVecRef<double, BandedMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<double, BandedMatrix>::VecType SxVecCompute<Simd>::divide (const double &y, const SxVecRef<double, BandedMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::subtract (const SxComplex16 &y, const SxVecRef<SxComplex16, Compact> &x) ;
SXVEC_INSTANCE SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::divide (const SxComplex16 &y, const SxVecRef<SxComplex16, Compact> &x) ;
SXVEC_INSTANCE SxVecResult<SxComplex16, PackedSymMatrix>::VecType SxVecCompute<Simd>::subtract (const SxComplex16 &y, const SxVecRef<SxComplex16, PackedSymMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<SxComplex16, PackedSymMatrix>::VecType SxVecCompute<Simd>::divide (const SxComplex16 &y, const SxVecRef<SxComplex16, PackedSymMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<SxComplex16, BandedMatrix>::VecType SxVecCompute<Simd>::subtract (const SxComplex16 &y, const SxVecRef<SxComplex16, BandedMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<SxComplex16, BandedMatrix>::VecType SxVecCompute<Simd>::divide (const SxComplex16 &y, const SxVecRef<SxComplex16, BandedMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<int, Compact>::VecType SxVecCompute<Simd>::subtract (const int &y, const SxVecRef<int, Compact> &x) ;
SXVEC_INSTANCE SxVecResult<int, Compact>::VecType SxVecCompute<Simd>::divide (const int &y, const SxVecRef<int, Compact> &x) ;
SXVEC_INSTANCE SxVecResult<int, PackedSymMatrix>::VecType SxVecCompute<Simd>::subtract (const int &y, const SxVecRef<int, PackedSymMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<int, PackedSymMatrix>::VecType SxVecCompute<Simd>::divide (const int &y, const SxVecRef<int, PackedSymMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<int, BandedMatrix>::VecType SxVecCompute<Simd>::subtract (const int &y, const SxVecRef<int, BandedMatrix> &x) ;
SXVEC_INSTANCE SxVecResult<int, BandedMatrix>::VecType SxVecCompute<Simd>::divide (const int &y, const SxVecRef<int, BandedMatrix> &x) ;
// --- Simd_scalar_subdiv
// --- The next 13 lines were generated from snippets/SxVecCompute_instantiate snippet UseIt_scalar_subdiv
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const double &y, const SxVecRef<double, SubMatrix>  &x) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const double &y, const SxVecRef<double, SubMatrix>  &x) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const double &y, const SxVecRef<double, Strided>  &x) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const double &y, const SxVecRef<double, Strided>  &x) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::subtract (const double &y, const SxVecRef<double, GeneralMatrix>  &x) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<UseIterator>::divide (const double &y, const SxVecRef<double, GeneralMatrix>  &x) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxComplex16 &y, const SxVecRef<SxComplex16, SubMatrix>  &x) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxComplex16 &y, const SxVecRef<SxComplex16, SubMatrix>  &x) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxComplex16 &y, const SxVecRef<SxComplex16, Strided>  &x) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxComplex16 &y, const SxVecRef<SxComplex16, Strided>  &x) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::subtract (const SxComplex16 &y, const SxVecRef<SxComplex16, GeneralMatrix>  &x) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<UseIterator>::divide (const SxComplex16 &y, const SxVecRef<SxComplex16, GeneralMatrix>  &x) ;
// --- UseIt_scalar_subdiv
// --- The next 41 lines were generated from snippets/SxVecCompute_instantiate snippet Simd_scalar
SXVEC_INSTANCE SxVecResult<float,Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE SxVecResult<float,Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE SxVecResult<float,Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE SxVecResult<float,Compact>::VecType SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<float, Compact> &x, const float &y);
SXVEC_INSTANCE SxVecResult<double,Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE SxVecResult<double,Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE SxVecResult<double,Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE SxVecResult<double,Compact>::VecType SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<double, Compact> &x, const double &y);
SXVEC_INSTANCE SxVecResult<SxComplex8,Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE SxVecResult<SxComplex8,Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE SxVecResult<SxComplex8,Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE SxVecResult<SxComplex8,Compact>::VecType SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<SxComplex8, Compact> &x, const SxComplex8 &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,Compact>::VecType SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<SxComplex16, Compact> &x, const SxComplex16 &y);
SXVEC_INSTANCE SxVecResult<int,Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<int, Compact> &x, const int &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<int, Compact> &x, const int &y);
SXVEC_INSTANCE SxVecResult<int,Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<int, Compact> &x, const int &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<int, Compact> &x, const int &y);
SXVEC_INSTANCE SxVecResult<int,Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<int, Compact> &x, const int &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<int, Compact> &x, const int &y);
SXVEC_INSTANCE SxVecResult<int,Compact>::VecType SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<int, Compact> &x, const int &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<int, Compact> &x, const int &y);
// --- Simd_scalar
// --- The next 17 lines were generated from snippets/SxVecCompute_instantiate snippet Simd_scalar_packed
SXVEC_INSTANCE SxVecResult<double,PackedSymMatrix>::VecType SxVecCompute<Simd>::add (const SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE SxVecResult<double,PackedSymMatrix>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE SxVecResult<double,PackedSymMatrix>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE SxVecResult<double,PackedSymMatrix>::VecType SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<double, PackedSymMatrix> &x, const double &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,PackedSymMatrix>::VecType SxVecCompute<Simd>::add (const SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,PackedSymMatrix>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,PackedSymMatrix>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
SXVEC_INSTANCE SxVecResult<SxComplex16,PackedSymMatrix>::VecType SxVecCompute<Simd>::Divide<true>::divide (const SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
SXVEC_INSTANCE void SxVecCompute<Simd>::Divide<true>::divideInPlace (SxVecRef<SxComplex16, PackedSymMatrix> &x, const SxComplex16 &y);
// --- Simd_scalar_packed
// --- The next 97 lines were generated from snippets/SxVecCompute_instantiate snippet SimdIt
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (const double* xPtr, SxVecConstIt<double, SubMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (SxVecConstIt<double, SubMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::add (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::add (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (const double* xPtr, SxVecConstIt<double, SubMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (SxVecConstIt<double, SubMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::subtract (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::subtract (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (const double* xPtr, SxVecConstIt<double, SubMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (SxVecConstIt<double, SubMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::multiply (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::multiply (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (const double* xPtr, SxVecConstIt<double, SubMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (SxVecConstIt<double, SubMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::divide (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::divide (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (const double* xPtr, SxVecConstIt<double, Strided> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (SxVecConstIt<double, Strided> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::add (const SxVecRef<double, Compact> &x, const SxVecRef<double, Strided>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::add (const SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (const double* xPtr, SxVecConstIt<double, Strided> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (SxVecConstIt<double, Strided> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::subtract (const SxVecRef<double, Compact> &x, const SxVecRef<double, Strided>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::subtract (const SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (const double* xPtr, SxVecConstIt<double, Strided> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (SxVecConstIt<double, Strided> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::multiply (const SxVecRef<double, Compact> &x, const SxVecRef<double, Strided>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::multiply (const SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (const double* xPtr, SxVecConstIt<double, Strided> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (SxVecConstIt<double, Strided> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::divide (const SxVecRef<double, Compact> &x, const SxVecRef<double, Strided>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::divide (const SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (const double* xPtr, SxVecConstIt<double, GeneralMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (SxVecConstIt<double, GeneralMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::add (const SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::add (const SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (const double* xPtr, SxVecConstIt<double, GeneralMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (SxVecConstIt<double, GeneralMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::subtract (const SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::subtract (const SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (const double* xPtr, SxVecConstIt<double, GeneralMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (SxVecConstIt<double, GeneralMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::multiply (const SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::multiply (const SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (const double* xPtr, SxVecConstIt<double, GeneralMatrix> it, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (SxVecConstIt<double, GeneralMatrix> it, const double* xPtr, ssize_t n, double* resPtr) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::divide (const SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdIt>::divide (const SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, SubMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (SxVecConstIt<SxComplex16, SubMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::add (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::add (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, SubMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (SxVecConstIt<SxComplex16, SubMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::subtract (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::subtract (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, SubMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (SxVecConstIt<SxComplex16, SubMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::multiply (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::multiply (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, SubMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (SxVecConstIt<SxComplex16, SubMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::divide (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::divide (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, Strided> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (SxVecConstIt<SxComplex16, Strided> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::add (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::add (const SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, Strided> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (SxVecConstIt<SxComplex16, Strided> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::subtract (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::subtract (const SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, Strided> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (SxVecConstIt<SxComplex16, Strided> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::multiply (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::multiply (const SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, Strided> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (SxVecConstIt<SxComplex16, Strided> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::divide (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::divide (const SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, GeneralMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::add (SxVecConstIt<SxComplex16, GeneralMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::add (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::add (const SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, GeneralMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::subtract (SxVecConstIt<SxComplex16, GeneralMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::subtract (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::subtract (const SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, GeneralMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::multiply (SxVecConstIt<SxComplex16, GeneralMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::multiply (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::multiply (const SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (const SxComplex16* xPtr, SxVecConstIt<SxComplex16, GeneralMatrix> it, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE void SxVecCompute<SimdIt>::divide (SxVecConstIt<SxComplex16, GeneralMatrix> it, const SxComplex16* xPtr, ssize_t n, SxComplex16* resPtr) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::divide (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix>  &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdIt>::divide (const SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
// --- SimdIt
// --- The next 25 lines were generated from snippets/SxVecCompute_instantiate snippet SimdCol
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::add (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::subtract (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::multiply (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::divide (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::add (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::subtract (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::multiply (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::divide (const SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::add (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::subtract (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::multiply (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<double> SxVecCompute<SimdCol>::divide (const SxVecRef<double, SubMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::add (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::subtract (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::multiply (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::divide (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::add (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::subtract (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::multiply (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::divide (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::add (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::subtract (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::multiply (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE SxVector<SxComplex16> SxVecCompute<SimdCol>::divide (const SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
// --- SimdCol
// --- The next 80 lines were generated from snippets/SxVecCompute_instantiate snippet Simd
SXVEC_INSTANCE typename SxVecResult<float, Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::add (float* resPtr, ssize_t, const float* x, const float *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<float, Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::subtract (float* resPtr, ssize_t, const float* x, const float *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<float, Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::multiply (float* resPtr, ssize_t, const float* x, const float *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<float, Compact>::VecType SxVecCompute<Simd>::divide (const SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::divide (float* resPtr, ssize_t, const float* x, const float *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::divideInPlace (SxVecRef<float, Compact> &x, const SxVecRef<float, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::add (double* resPtr, ssize_t, const double* x, const double *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::subtract (double* resPtr, ssize_t, const double* x, const double *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::multiply (double* resPtr, ssize_t, const double* x, const double *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::divide (const SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::divide (double* resPtr, ssize_t, const double* x, const double *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::divideInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex8, Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::add (SxComplex8* resPtr, ssize_t, const SxComplex8* x, const SxComplex8 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex8, Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::subtract (SxComplex8* resPtr, ssize_t, const SxComplex8* x, const SxComplex8 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex8, Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::multiply (SxComplex8* resPtr, ssize_t, const SxComplex8* x, const SxComplex8 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex8, Compact>::VecType SxVecCompute<Simd>::divide (const SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::divide (SxComplex8* resPtr, ssize_t, const SxComplex8* x, const SxComplex8 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::divideInPlace (SxVecRef<SxComplex8, Compact> &x, const SxVecRef<SxComplex8, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::add (SxComplex16* resPtr, ssize_t, const SxComplex16* x, const SxComplex16 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::subtract (SxComplex16* resPtr, ssize_t, const SxComplex16* x, const SxComplex16 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::multiply (SxComplex16* resPtr, ssize_t, const SxComplex16* x, const SxComplex16 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::divide (const SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::divide (SxComplex16* resPtr, ssize_t, const SxComplex16* x, const SxComplex16 *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::divideInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<int, Compact>::VecType SxVecCompute<Simd>::add (const SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::add (int* resPtr, ssize_t, const int* x, const int *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::addInPlace (SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<int, Compact>::VecType SxVecCompute<Simd>::subtract (const SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::subtract (int* resPtr, ssize_t, const int* x, const int *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::subtractInPlace (SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<int, Compact>::VecType SxVecCompute<Simd>::multiply (const SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::multiply (int* resPtr, ssize_t, const int* x, const int *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::multiplyInPlace (SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;

SXVEC_INSTANCE typename SxVecResult<int, Compact>::VecType SxVecCompute<Simd>::divide (const SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::divide (int* resPtr, ssize_t, const int* x, const int *y);
SXVEC_INSTANCE void SxVecCompute<Simd>::divideInPlace (SxVecRef<int, Compact> &x, const SxVecRef<int, Compact> &y) ;
// --- Simd
// --- The next 9 lines were generated from snippets/SxVecCompute_instantiate snippet UseIt_InPlace_mixed
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Compact>  &x,
      const SxVecRef<double, Compact> &y);
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Compact>  &x,
      const SxVecRef<double, Compact> &y);
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Compact>  &x,
      const SxVecRef<double, Compact> &y);
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, Compact>  &x,
      const SxVecRef<double, Compact> &y);
// --- UseIt_InPlace_mixed
// --- The next 97 lines were generated from snippets/SxVecCompute_instantiate snippet UseIt_inPlace
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, Compact> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, SubMatrix> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, Strided> &x, const SxVecRef<double, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<double, GeneralMatrix> &x, const SxVecRef<double, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, Compact> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, SubMatrix> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, Strided> &x, const SxVecRef<SxComplex16, GeneralMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Compact> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, SubMatrix> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::addInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::subtractInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::multiplyInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
SXVEC_INSTANCE void SxVecCompute<UseIterator>::divideInPlace (SxVecRef<SxComplex16, GeneralMatrix> &x, const SxVecRef<SxComplex16, Strided> &y) ;
// --- UseIt_inPlace
// --- The next 9 lines were generated from snippets/SxVecCompute_instantiate snippet Simd_func
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::sqr (const SxVecRef<double, Compact> &vec) ;
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::sqrInPlace (SxVecRef<double, Compact> &&vec) ;
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::cub (const SxVecRef<double, Compact> &vec) ;
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::cubInPlace (SxVecRef<double, Compact> &&vec) ;
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::exp (const SxVecRef<double, Compact> &vec) ;
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::expInPlace (SxVecRef<double, Compact> &&vec) ;
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::sqrt (const SxVecRef<double, Compact> &vec) ;
SXVEC_INSTANCE typename SxVecResult<double, Compact>::VecType SxVecCompute<Simd>::sqrtInPlace (SxVecRef<double, Compact> &&vec) ;
// --- Simd_func
// --- The next 3 lines were generated from snippets/SxVecCompute_instantiate snippet Simd_transpose
SXVEC_INSTANCE void SxVecCompute<Simd>::transpose (ssize_t N, ssize_t M, const double* src, ssize_t strideS, double* dest, ssize_t strideD) ;
SXVEC_INSTANCE void SxVecCompute<Simd>::transpose (ssize_t N, ssize_t M, const SxComplex16* src, ssize_t strideS, SxComplex16* dest, ssize_t strideD) ;
// --- Simd_transpose
// --- The next 3 lines were generated from snippets/SxVecCompute_instantiate snippet Simd_conj
SXVEC_INSTANCE typename SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::conj (const SxVecRef<SxComplex16, Compact> &vec) ;
SXVEC_INSTANCE typename SxVecResult<SxComplex16, Compact>::VecType SxVecCompute<Simd>::conjInPlace (SxVecRef<SxComplex16, Compact> &&vec) ;
// --- Simd_conj
// --- The next 2 lines were generated from snippets/SxVecCompute_instantiate snippet Simd_adjoint
SXVEC_INSTANCE void SxVecCompute<Simd>::adjoint (ssize_t N, ssize_t M, const SxComplex16* src, ssize_t strideS, SxComplex16* dest, ssize_t strideD) ;
// --- Simd_adjoint
