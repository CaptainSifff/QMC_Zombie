#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

#ifndef GCC_VER
#define GCC_VER(x,y,z)  ((x) * 10000 + (y) * 100 + (z))
#endif

#include <cmath>
/**
Linear Interpolation between l and h. a must be from the interval [0,1]
@param a how much to interpolate
@param l the lower point
@param h the higher point
*/
template <typename T, typename U>
inline T lerp(U a , T l, T h) throw()
{
    return l + (h-l) * a;
}

template <typename FPType, bool FT_is_real = false>
struct MatsubaraFT
{
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (GCC_VERSION > GCC_VER(4,5,0))
    static constexpr FPType normalizationfactor = 1.0;
#else
    static const FPType normalizationfactor = 1.0;
#endif
    template<class GF, typename OutputType>
    static inline void call(GF& gf, FPType omegam, const FPType a, const FPType b, OutputType& afac, OutputType& bfac)
    {
      FPType invomegam = 1.0/omegam;
      typedef typename GF::RetType Type;
      Type gomegaplus = gf(omegam);
      Type gomegaminus = gf(-omegam);
      afac = gomegaplus + gomegaminus - 2.0 * b * invomegam*invomegam;
      bfac = Type(0.0, -1.0)*(gomegaplus - gomegaminus - Type(0.0, 2.0*a*invomegam));
    }
};

template <typename FPType>
struct MatsubaraFT<FPType, true>
{
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (GCC_VERSION > GCC_VER(4,5,0))
    static constexpr FPType normalizationfactor = 2.0;
#else
    static const FPType normalizationfactor = 2.0;
#endif
    template<class GF, typename OutputType>
    static inline void call(GF& gf, FPType omegam, const FPType a, const FPType b, OutputType& afac, OutputType& bfac)
    {
      FPType invomegam = 1.0/omegam;
      typename GF::RetType gp = gf(omegam);
      bfac = imag(gp) - a*invomegam;
      afac = real(gp) - b*invomegam*invomegam;
    }
};

template <class GF, typename OutputType, typename FPType>
inline void matsubarafouriertransform(GF& gf, OutputType* g, FPType beta, FPType betaslice, unsigned int glen, unsigned int betaslices)
{
  /**
  Do not try to parallelize anything here using OpenMP
  both for loops can then potentially lead to race conditions.
  The only possibility would be to synchronize the access to g[k] in the update relation using an atomic operation.
  OpenMP does not support atomic operations for overloaded operators.
  And the += operator is overloaded for the complex<> type.
  A possible way out would be the use of the C++11 thread library. 
  But only gcc-4.7 can report the maximum number of hardware supported threads. 
  */
    for (unsigned int k = 0; k < glen; ++k) g[k] = 0;
    const unsigned int m_max = 30000;//Nr of Matsubara frequencies.
    //We implement the Large frequency behaviour trick of David/Fakher
    const unsigned int nr_interesting_points = 5;//obvious
    FPType a = 0.0;
    FPType b = 0.0;
    typedef typename GF::RetType Type;
    for ( int k = m_max; k > static_cast<int>(m_max - nr_interesting_points); --k)
    {
        FPType omegam = static_cast<FPType>(M_PIl) * static_cast<FPType>(2*k+1)/beta;
        Type gom = gf(omegam);
        b += real(gom) * omegam * omegam;
        a += imag(gom) * omegam;
    }
    b /= static_cast<FPType>(nr_interesting_points);
    a /= static_cast<FPType>(nr_interesting_points);
    for (uint m = 0; m <= m_max; ++m)
    {
        FPType com = static_cast<FPType>(M_PIl * static_cast<FPType>(2*m+1));
        FPType omegam = com/beta;
	OutputType afac, bfac;
	MatsubaraFT<FPType, GF::has_real_FourierTransform>::call(gf, omegam, a, b, afac, bfac );
        FPType fac = com/betaslices;
        FPType temp = sin(fac/2.0);
        const FPType alpha = 2.0*temp*temp;
        const FPType betafac = sin(fac);
        FPType cn = 1.0;
        FPType sn = 0.0;
        for (unsigned int k = 0; k < glen; ++k)
        {
            g[k] += cn * afac + sn * bfac;
            //update the recursion relation for the trigonometric functions
            const FPType tcn = cn;
            cn -= (alpha * cn + betafac * sn);
            sn -= (alpha * sn - betafac * tcn);
        }
    }
    for(unsigned int k = 0; k < glen; ++k)
    {
        g[k] *= MatsubaraFT<FPType, GF::has_real_FourierTransform>::normalizationfactor/beta;
        g[k] += 0.5;
        g[k] -= b/2.0*(k * betaslice - beta/2.0); 
    }
}
#endif