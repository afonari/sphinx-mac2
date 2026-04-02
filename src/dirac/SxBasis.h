// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institut fuer Eisenforschung GmbH
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxrepo.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#ifndef _SX_BASIS_H_
#define _SX_BASIS_H_

#include <SxMemConsumer.h>
#include <SxDiracLib.h>
#include <SxString.h>
#include <SxArray.h>
#include <SxPrecision.h>
#include <SxVector.h>

/* --- First-class basis sets
   These basis sets are directly supported by SxBasis. Additional basis sets
   are handled via the SxExtraBasis<> interface.
*/
class SxGBasis;
class SxRBasis;
class SxAOBasis;
class SxRadBasis;
class SxRadialBasis;

/*  projectTo (B, v) projects vector v from its current basis to basis B;
                     the result is a reference. This allows to do cheap
                     "projections" to the same basis by returning a reference
                     to the current data.
                     This is the function behind (B | v)
    toBasis (B,v)    also projects vector v to basis B, but the return type is
                     a submatrix reference. This allows to extract a reference
                     to part of a composite basis
                     (e.g. PAWBasis = G basis + partial wave basis)
                     This is the function behind   v.to (B) .

                     This function is not overloaded by REGISTER_PROJECTOR, but
                     must be explicitly overloaded.

                     NOTE: at present, this is restricted to double and
                     SxComplex16 because we do not use single-precision
                     composite basis sets anywhere.
   
    note: RESTYPE must be the same as B::CoeffType, but we cannot
    use that here because the SxBasis base class must be defined
    before any of the derived basis classes
*/
#define VIRTUAL_PROJECT_TO(B,RESTYPE)                                   \
      typedef SxVecRef<RESTYPE> SxVecRef_ ## B;                         \
      virtual SxVecRef<RESTYPE>                                         \
      projectTo (const B *, const SxVecRef<float> &) const              \
      {                                                                 \
         projectionFailed ( #B , "(float->" #RESTYPE ")");              \
         return SxVecRef<RESTYPE>();                                    \
      }                                                                 \
      virtual SxVecRef<RESTYPE>                                         \
      projectTo (const B *, const SxVecRef<double> &) const             \
      {                                                                 \
         projectionFailed ( #B , "(double->" #RESTYPE ")");             \
         return SxVecRef<RESTYPE>();                                    \
      }                                                                 \
      virtual SxVecRef<RESTYPE>                                         \
      projectTo (const B *, const SxVecRef<SxComplex8> &) const         \
      {                                                                 \
         projectionFailed ( #B , "(SxComplex8->" #RESTYPE ")");         \
         return SxVecRef<RESTYPE>();                                    \
      }                                                                 \
      virtual SxVecRef<RESTYPE>                                         \
      projectTo (const B *, const SxVecRef<SxComplex16> &) const        \
      {                                                                 \
         projectionFailed ( #B , "(SxComplex16->" #RESTYPE ")");        \
         return SxVecRef<RESTYPE>();                                    \
      }                                                                 \
                                                                        \
      virtual SxVecRef<RESTYPE, SubMatrix>                              \
      toBasis (const B *b, const SxVecRef<double> &v) const             \
      {                                                                 \
         SxVecRef<RESTYPE> res = projectTo (b, v);                      \
         SxVecRef<RESTYPE,SubMatrix> resSub = res.getRef<SubMatrix> (0, \
               res.getNRows (), res.getNCols (), 1, res.getNRows ());   \
         resSub.auxData = res.auxData;                                  \
         return resSub;                                                 \
      }                                                                 \
                                                                        \
      virtual SxVecRef<RESTYPE, SubMatrix>                              \
      toBasis (const B *b, const SxVecRef<SxComplex16> &v) const        \
      {                                                                 \
         SxVecRef<RESTYPE> res = projectTo (b, v);                      \
         SxVecRef<RESTYPE,SubMatrix> resSub = res.getRef<SubMatrix> (0, \
               res.getNRows (), res.getNCols (), 1, res.getNRows ());   \
         resSub.auxData = res.auxData;                                  \
         return resSub;                                                 \
      }

/* This macro implements the virtual projectTo interface for a projector
   to basis "TO". Make sure that the target basis has typedef'ed
   CoeffType (the vector element (scalar) type) in agreement with the
   result of PROJECTOR and with the VIRTUAL_PROJECT_TO list inside
   the SxBasis definition below.
   */
#define REGISTER_PROJECTOR(TO,PROJECTOR)                                \
     virtual SxVecRef_ ## TO                                            \
     projectTo (const TO *b, const SxVecRef<float> &in) const           \
     {                                                                  \
        return PROJECTOR (b, in);                                       \
     }                                                                  \
     virtual SxVecRef_ ## TO                                            \
     projectTo (const TO *b, const SxVecRef<double> &in) const          \
     {                                                                  \
        return PROJECTOR (b, in);                                       \
     }                                                                  \
     virtual SxVecRef_ ## TO                                            \
     projectTo (const TO *b, const SxVecRef<SxComplex8> &in) const      \
     {                                                                  \
        return PROJECTOR (b, in);                                       \
     }                                                                  \
     virtual SxVecRef_ ## TO                                            \
     projectTo (const TO *b, const SxVecRef<SxComplex16> &in) const     \
     {                                                                  \
        return PROJECTOR (b, in);                                       \
     }

#define FOR_ALL_VECTYPES(MACRO) \
   MACRO(float) \
   MACRO(double) \
   MACRO(SxComplex8) \
   MACRO(SxComplex16)

#define VIRTUAL_PROJECT_FROM(BASIS, COEFF_TYPE)                         \
virtual SxVecRef<CoeffType> projectFrom (const BASIS *fromBasis,        \
                                 const SxVecRef<COEFF_TYPE> &) const    \
{                                                                       \
   SX_CHECK (fromBasis);                                                \
   projectionFailed (#BASIS,dynamic_cast<const SxBasis*>(this));        \
   return SxVecRef<CoeffType> ();                                       \
}

#define PROJECT_TOEXTRA_VIA_FROM(TARGET, TVEC)             \
virtual SxVecRef<TARGET::CoeffType>                        \
projectTo (const SxExtraBasis<TARGET::CoeffType> *to,      \
           const SxVecRef<TVEC> &vec) const                \
{                                                          \
   const TARGET *target = dynamic_cast<const TARGET*>(to); \
   SX_CHECK (target);                                      \
   return target->projectFrom (this, vec);                 \
}

template <class T> class SxExtraBasis;

/** \brief Abstract basis-set class used for direct Dirac-like projections.

    \b SxBasis = S/PHI/nX Basis-sets 

    This class is a major class to establish the \ref page_dirac. It is 
    used to perform direct projections onto basis-sets. Direct means that
    the destination basis is explicitly known and Dirac's bra-ket notation
    is used.

    \ingroup group_dirac
    \sa      \ref page_dirac
    \sa      SxBasis
    \author  Sixten Boeck
  */
class SX_EXPORT_DIRAC SxBasis
: public SxMemConsumer,
  public SxThis<SxBasis>
{
   private:
      /** \brief Pointers to registered bases

          It is sometimes necessary that basis transformations
          require auxiliary data. These are stored in one of the
          bases (and should be mutable). In order to relate this
          data correctly with the other basis, a pointer to the
          second basis is registered and stored here.

          To add a pointer use the function ::registerBasis.

          If several bases of the same type are registered, the
          id of a specific basis can be obtained via getBasisId.

          \sa \ref page_dirac */
      mutable SxArray<const SxBasis *> registeredBases;

   public:
      /// \brief Number of sampling points of the corresponding Dirac vector.
      virtual ssize_t getNElements () const = 0;

      /// \brief standard destructor
      virtual ~SxBasis ();

      /** \brief Register basis unless already registered
          @param basisPtr basis to be registered
          @return true if basis was new and registered,
                  false if basis was already registered
          @note This function needs rarely to be called outside the basis
                layer, since bases may autoregister in the first projection.
                If newly implemented bases require registration at some point,
                they should autoregister right in place. (That's why we made
                this stuff mutable).
        */
      bool registerBasis (const SxBasis &basis) const;


      /** \brief Register a basis of unknown type

          \note This is a hook-in for derived bases that
                provide special-case registration. The correct type
                must be obtained by dynamic_cast.
          \note The deregistration hook-in is deregister.
        */
      virtual void registerUnknown (const SxBasis &basis) const
      {
         registerBasis (basis);
      }

   protected:
      /** \brief Deregister a basis
        This performs the SxBasis deregistration process.
       */
      void deregisterBasis (const SxBasis *) const;
      
      /** \brief Deregister a basis from derived basis
        This function needs to be overloaded if the deregistration
        should remove some internal data of the derived basis.
        Use dynamic_cast<> to figure out the type of the deregistered 
        basis.

        Remember to call deregisterBasis at the end if you overload
        this.
          
       */
      virtual void deregister (const SxBasis *basis) const
      {
         deregisterBasis (basis);
      }
      
      /** Deregister all bases
          \note This must be called in the derived basis' destructor.
        */
      void deregisterAll ();
      
   public:

      /** \brief Get basis id
          If more than one basis of some type is registered, each
          basis gets an id, starting from 0. In this way, member data
          of the derived basis can be correctly be mapped to the basis
          in question.
          @return the id, or -1 if basis is not registered

          \example
          \code
int iFFT = getBasisId<SxRBasis> (targetBasis);
fft3d(iFFT).fftForward (data, res);
          \endcode
        */
      template <class BasisType>
      int getBasisId (const BasisType* basis) const;

      /// Check if a basis is registered
      bool isRegistered (const SxBasis *basis) const
      {
         for (int i = 0; i < registeredBases.getSize (); ++i)
            if (registeredBases(i) == basis) return true;
         return false;
      }

      /** \brief placeholder for any \f$
                    \mathrm{tr}_\mathbf{B}(\hat{\varrho}\hat{\mathbf{A}})
                 \f$  */
#define VIRTUAL_TR(T)                                                        \
      virtual Real8 tr (const SxVecRef<T> &) const                           \
      {                                                                      \
         printf ("Trace for SxVecRef<%s> in %s basis not yet implemented.\n",\
                 #T, getType ().ascii ());                                   \
         SX_EXIT; return 0.;                                                 \
      }
      FOR_ALL_VECTYPES (VIRTUAL_TR)

   
   private:
      /// Crash on an unimplemented projection
      void projectionFailed (const char *basis, const char *types) const;
   public:
      /// Crash on an unimplemented projection
      static void projectionFailed (const char* from, const SxBasis* to);

      /** \brief placeholder for any <??|X> projector */
#define VIRTUAL_TOANY(TYPE)                                        \
      virtual SxVecRef<TYPE>                                       \
      projectTo (const SxBasis *, const SxVecRef<TYPE> &) const    \
      {                                                            \
         projectionFailed ( "SxBasis" , "(" #TYPE "->" #TYPE ")"); \
         return SxVecRef<TYPE> ();                                 \
      }
      FOR_ALL_VECTYPES (VIRTUAL_TOANY)

      /** \brief placeholder for any <R|X> projector */
      VIRTUAL_PROJECT_TO (SxRBasis, SxComplex16)
      /** \brief placeholder for any <G|X> projector */
      VIRTUAL_PROJECT_TO (SxGBasis, SxComplex16)
      /** \brief placeholder for any <mu|X> projector */
      VIRTUAL_PROJECT_TO (SxAOBasis, SxComplex16)
      /** \brief placeholder for any <r|X> projector */
      VIRTUAL_PROJECT_TO (SxRadBasis, double)
      /** \brief placeholder for any <radial|X> projector */
      VIRTUAL_PROJECT_TO (SxRadialBasis, double)

#define VIRTUAL_PROJECT_TO_EXTRA(TBAS,TVEC)                              \
      virtual SxVecRef<TBAS> projectTo (const SxExtraBasis<TBAS> *,      \
                                        const SxVecRef<TVEC> &) const    \
      {                                                                  \
         projectionFailed ("SxExtraBasis<" #TBAS ">", #TVEC "->" #TBAS); \
         return SxVecRef<TBAS> ();                                       \
      }
      VIRTUAL_PROJECT_TO_EXTRA (double,float)
      VIRTUAL_PROJECT_TO_EXTRA (double,double)
      VIRTUAL_PROJECT_TO_EXTRA (SxComplex16,float)
      VIRTUAL_PROJECT_TO_EXTRA (SxComplex16,double)
      VIRTUAL_PROJECT_TO_EXTRA (SxComplex16,SxComplex8)
      VIRTUAL_PROJECT_TO_EXTRA (SxComplex16,SxComplex16)

      /** \brief Print debug information about the basis */
      virtual void print () const;

      /** Very simple description of basis a la "|R>" or "|G>"
        */
      virtual SxString getType () const = 0;

      /// Scalar product
      virtual double scalarProduct (const SxVecRef<double> &x,
                                    const SxVecRef<double> &y) const;

      /// Scalar product
      virtual SxComplex16 scalarProduct (const SxVecRef<SxComplex16> &x,
                                         const SxVecRef<SxComplex16> &y) const;

   protected:

      /** \brief Standard constructor */
      SxBasis () {/* empty */}

   public:

      
};

// standalone for getNElements to allow SxVector (const SxBasis &);
/// \brief Number of sampling points of the corresponding Dirac vector.
inline ssize_t getNElements (const SxBasis &b)
{
   return b.getNElements ();
}

template <class BasisType>
int SxBasis::getBasisId (const BasisType* basis) const
{
   int i, id = 0;
   // find basis of desired type
   for (i = 0; i < registeredBases.getSize (); ++i)  {
      if (registeredBases(i) == basis) return id;
      // increment id, if other basis is found
      if (dynamic_cast<const BasisType*>(registeredBases(i))) id++;
   }
   // basis is not registered
   SX_EXIT;
}

/** < B | vec >
  \sa     \ref page_dirac
  \author Sixten Boeck
  */
template<class T, class B>
SxVecRef<typename B::CoeffType>
operator| (const B &b, const SxVecRef<T> &v)
{
   SX_CHECK (static_cast<const SxBasis*> (&b)); // must be a basis
   SX_CHECK (v.getBasisPtr());
   return v.getBasisPtr()->projectTo (&b, v);
}

/** < B(?) | vec >
  \sa     \ref page_dirac
  \author C. Freysoldt
  @note: This is the projection operator for an anonymous basis (type not
         known at compile time). The return type is the same as the
         input type! Don't use if this is not appropriate.
  */
template<class T>
SxVecRef<T>
operator| (const SxBasis &b, const SxVecRef<T> &v)
{
   SX_CHECK (v.getBasisPtr());
   SxVecRef<T> res;
   return v.getBasisPtr()->projectTo (&b, v);
}

/** \brief Basis class for additional basis sets

    Additional basis sets outside dirac/ must be derived from SxExtraBasis<T>.
    The template parameter defines the "standard" coefficient type that results
    when projecting to the basis.

    If a group of mutually convertible basis sets is created, one should
    introduce a common base class that provides additional VIRTUAL_PROJECT_FROM
    slots. The transformations should then be interfaced in the target basis
    via the new projectFrom interfaces.
    Then use the PROJECT_TOEXTRA_VIA_FROM macro in the derived classes to create the
    basis-type resolving indirection
    (SxBasis*)Y->projectTo ((SxExtraBasis<T>*)X, vec(Y)) -> X->projectFrom (Y, vec(Y)) .

  */
template <class T> class SxExtraBasis;
template <>
class SxExtraBasis<SxComplex16> : public SxBasis {
   public:
      typedef SxComplex16 CoeffType;
      VIRTUAL_PROJECT_FROM (SxGBasis, SxComplex16);
      //VIRTUAL_PROJECT_FROM (SxRBasis, SxComplex16); // not needed
      VIRTUAL_PROJECT_FROM (SxAOBasis, SxComplex16);
};
template <>
class SxExtraBasis<double> : public SxBasis {
   public:
      typedef double CoeffType;
      VIRTUAL_PROJECT_FROM (SxRadialBasis, double);
};

// --------------------------------------------------------------------------
/** \brief placeholder for \f$ \mathrm{tr} (\hat{\varrho} \hat{A}) \f$

    Dependent on the basis-set this function returns the expectation
    value of an operator \f$ \langle \hat{A} \rangle \f$.
    Using the notation of a density matrix the expectation value can
    be evaluated according to
    \f[
       \langle \hat{A} \rangle
         =
       \mathrm{tr} (\hat{\varrho} \hat{A})
    \f]
    with the density matrix
    \f[
       \hat{\varrho} = \sum_i | \Psi_i \rangle  \langle \Psi_i |
    \f]

    The actual implementation does not use a density matrix. Instead the
    product of the density with an operator is taken as argument which
    is a Dirac vector again. This Dirac vector is aware of the basis-set
    and calls the corresponding ::tr function of this basis-set.
    Each basis-set should overload the corresponding ::tr function.
    Since the ::tr function is basis-set dependent we write the
    basis-set as index
    \f[
       \langle \hat{A_{\mathbf{B}}} \rangle
         =
       \mathrm{tr}_{\mathbf{B}} (\hat{\varrho} \hat{A})
    \f]
    with \b B beeing the basis-set.
    \param  in  the evaluated product of \f$ \hat{\varrho} \hat{A} \f$
    \return     \f$ \mathrm{tr} (\hat{\varrho} \hat{A} ) \f$*/
template<class T>
Real8 tr (const SxVecRef<T> &in)
{
   SX_CHECK (in.getBasisPtr());
   return in.getBasisPtr()->tr (in);
}

template <class T>
T
operator| (const SxVecRef<T> &x, const SxVecRef<T> &y)
{
   return x.getBasisPtr ()->scalarProduct (x, y);
}


#endif /* _SX_BASIS_H_ */
