/***************************************************************************
 *   Copyright (C) 2014 by Florian Goth   *
 *   fgoth@wthp095   *
 *                                                                         *
 *   All rights reserved.                                                  *
 *                                                                         *
 *   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: *
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. *
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. *
 *     * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. *
 *                                                                         *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   *
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     *
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR *
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR *
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, *
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,   *
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    *
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF *
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING  *
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 *   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          *
 ***************************************************************************/
#ifndef POWERLAW_RASHBA_CHAIN_H
#define POWERLAW_RASHBA_CHAIN_H
#include <complex>
#include <limits>
#include <tr1/cmath>
#include "ddqmc.h"
#include "Vertex.h"
#include "Greensfunction.h"
/**
An implementation of a 1D chain with an added Rashba-type interaction.
the transition rates decay power-law like
I assume that the lattice constant is a = 1
Due to time-reversal symmetry and the particular form of the coefficients of the Hamiltonian,
It always returns real values.
*/
template<typename FPType = float>
class PowerLawRashbaChain
{
  public:
    enum {timeevolution = 0,
    has_Spin_symmetry = true,//helical base
    has_Giomega = false,
    Is_Impurity_model = false,
    has_TRS = true,
    };
    typedef Hubbard_Vertex<FPType> Vertex;
    typedef std::complex<FPType> FreeGreensFunctionReturnValueType;///< a typedef of the data type of the free greensfunction
    /**
    This evaluates the value of the free particle Green-function for the two given vertices
    @param v1 the first vertex
    @param v2 the second vertex
    @return the value of the free greensfunction evaluated with the given vertices
    */
    static inline FreeGreensFunctionReturnValueType eval(const Vertex& v1, const Vertex& v2) throw();
    /**
    A function for initializing the tables that make up the Greensfunction. Beta and the nr of sites are read from config files
    @param CFG a class that contains the necessary information to extract the parameters.
    */
    template <class CFG>
    static inline void init(CFG&);
    /**
    This function frees the memory used by the tables for the Greensfunction
    */
    static inline void tidyup();
    /**
    To access the nr of atoms of the chain
    @return the nr of atoms in the chain
    */
    static unsigned int getLen() throw()
    {
        return N;
    }
    /**
    Get the length of the interactioninterval.
    @return the inverse temperature beta
    */
    static FPType getContourLen() throw()
    {
        return beta;
    }
private:
    static std::complex<FPType>* g;///< the free particle greensfunction. We store the up and down part next to each other in a complex number. If we use the symmetry some time in the future this might be helpful. Note that we switch to a flat array.
    static unsigned int slices;///< the number of timeslices. This means the resolution of the free greensfunction in tau space
    static FPType betaslice;///< The length of one timeslice on the tau axis
    static FPType beta;///< The inverse temperature beta
    static unsigned int N;///< The length of the chain
    /**
    the Dispersionrelation -2t cos(k) -my . Note that the chemical potential is included here.
    @param k the momentum
    @param t the kinetic energy
    @param my the chemical potential
    @param lambda the spin-orbit strength
    @return the energy
    */
    static inline FPType disp(int k, FPType t, FPType my, FPType lambda) throw();
//     /**
//     * The free particle Greens function in k-space for the up - up part and for the down-down part
//     @param k the momentum
//     @param tau the imaginary time
//     @param t the kinetic energy
//     @param my the chemical potential
//     @param lambda the spin-orbit strength
//     */
//     static inline FPType GreenNullKSpacediag(FPType k, FPType tau, FPType t, FPType my, FPType lambda) throw();
//     /**
//     * The free particle Greens function in k-space for the up - down part as well as the down-up part
//     * since they are symmetric.
//     @param k the momentum
//     @param tau the imaginary time
//     @param t the kinetic energy
//     @param my the chemical potential
//     @param lambda the spin-orbit strength
//     */
//     static inline FPType GreenNullKSpaceoffdiag(FPType k, FPType tau, FPType t, FPType my, FPType lambda) throw();
    /**
     * This Green's function should suffice to specify the necessary Green's function.
     * After the Fourier-transform the spin-diagonal G can be recovered from the real part and the 
     * off-diagonal from the imaginary part.
    @param k the momentum
    @param tau the imaginary time
    @param t the kinetic energy
    @param my the chemical potential
    @param lambda the spin-orbit strength
    */
    static inline FPType GreenNullKSpace(int k, FPType tau, FPType t, FPType my, FPType lambda) throw();
    /**
     * This implements Li_s( exp(w)).
     * From here:  http://www.cs.kent.ac.uk/pubs/1992/110/
     * Note that there's a log in this series expansion. The imaginary therefore depends on the branch cut of the log.
     * 
     * FIXME: something is still wrong in here. but it's the same if the expression is evaluated with Mathematica
     * */
    static inline std::complex<FPType> PolyLog_Exp(const unsigned int s, std::complex<FPType> w)
    {//s>=2
      std::cout<<"Integer Series"<<std::endl;
      std::complex<FPType> res = std::tr1::__detail::__riemann_zeta(static_cast<FPType>(s));//optimization possibility: s are positive integers
      std::complex<FPType> wpower = w;
      FPType fac = 1.0;
      FPType harmonicN = 1.0;//HarmonicNumber_1
      for (uint k = 1; k <= s-2; ++k)
      {
	res += wpower*fac*std::tr1::__detail::__riemann_zeta(static_cast<FPType>(s - k));
	wpower *= w;
	FPType temp = 1.0/(1.0 + k);
	fac *= temp;
	harmonicN += temp;
//	std::cout<<k<<" "<<res<<std::endl;
      }
      //harmonicN now contains H_{s-1}
      //fac should be 1/(n-1)!
      res += (harmonicN - std::log(-w))*wpower*fac;
      wpower *= w;
      fac /= s;
      res -= wpower*fac/2.0;
      wpower *= w;
      //now comes the remainder of the series.
      const FPType tp = 2.0 * M_PI;
      const std::complex<FPType> pref = wpower/M_PI/tp;
      const unsigned int maxit = 200;
      unsigned int j = 1;
      bool terminate = false;
      fac /= (s+1.0);//(1/(n+1)!)
      //subtract the zeroth order term.
      res -= M_PI*M_PI/6.0*fac * pref;
//      std::cout<<res<<std::endl;
      //remainder of series
      fac *= 3.0*2.0/(s + 2.0)/(s+3.0);
      std::complex<FPType> upfac = -(w/tp)*(w/tp);
      std::complex<FPType> w2 = upfac;
      while (!terminate)//assume uniform convergence
      {
        FPType rzarg = static_cast<FPType>(2*j+2);
        FPType rz = std::tr1::__detail::__riemann_zeta(rzarg);
//        std::cout<<rz<<" "<<fac<<" "<<w2<<std::endl;
	std::complex<FPType> nextterm = (rz*fac)*w2;
	w2 *= upfac;
	fac *= rzarg/(rzarg + s) * (rzarg+1.0)/(rzarg + s + 1.0);
	++j;
	terminate = (fpequal( std::abs(res - pref*nextterm), std::abs(res) ) || (j > maxit));
	res -= pref * nextterm;
      }
   std::cout<<"Iterations in Integer Series: "<<j<<std::endl;
//       std::cout<<res2<<" "<<res2*pref<<std::endl;
      return res;
    }
    static inline std::complex<FPType> PolyLog_Exp_pos(const FPType s, std::complex<FPType> w)
    {
      std::complex<FPType> res =  std::tr1::__detail::__riemann_zeta(s);
      std::complex<FPType> wpower = w;
      FPType fac = 1.0;
      const int m = static_cast<int>(std::floor(s));
      for (uint k = 1; k <= m; ++k)
      {
	res += wpower*fac*std::tr1::__detail::__riemann_zeta(static_cast<FPType>(s - k));
	wpower *= w;
	FPType temp = 1.0/(1.0 + k);
	fac *= temp;
      }
      //fac should be 1/(m+1)!
      res += std::tgamma(1-s)*std::pow(-w, s-1);
      const FPType tp = 2.0 * M_PI;
      const FPType pref = 2.0 * std::pow(tp, s-1);
      //now comes the remainder of the series
      const unsigned int maxit = 100;
      unsigned int j = 0;
      bool terminate = false;
      std::complex<FPType> wup = w/tp;
      std::complex<FPType> w2 = std::pow(wup, m+1);
      std::complex<FPType> gam = std::tgamma(2.0-s+m)*fac; //here we factor up the ratio of Gamma(1 - s + k)/k! . This ratio should be well behaved even for large k
      FPType sp, cp;
      sincos(M_PI/2.0 * s, &sp, &cp);
      while (!terminate)//assume uniform convergence
      {//FIXME: optimize.
	int idx = m + 1 + j;
	FPType zetaarg = 1 - s + idx;
	FPType sine;
	if(idx & 1)//save the reperated calculation of the sines
	{/*odd*/
	  sine = cp;
	  if ( !((idx-1)/ 2 & 1) )
	    sine = -sine;
	}
	else
	{/*even*/
	  sine = sp;
	  if((idx/2) & 1)
	    sine = -sine;
	}
	std::complex<FPType> nextterm = (std::tr1::__detail::__riemann_zeta(zetaarg) * sine * gam) * w2;
//	std::cout<<j<<" "<<nextterm<<" used Gamma = "<<gam<<std::endl;
	w2 *= wup;
	gam *= zetaarg/(1.0 + idx);
	++j;
	terminate = (fpequal( std::abs(res + pref*nextterm), std::abs(res) ) || (j > maxit));
	res += pref * nextterm;
      }
      std::cout<<"Iterations in PolyLogExp_pos: "<<j<<std::endl;
      return res;
    }
    static inline std::complex<FPType> PolyLog_Exp_neg(const FPType s, std::complex<FPType> w)
    {//basic general loop, but s is a negative quantity here
      //TODO: optimize/fix for the case of negative Integers
      std::cout<<"general loop"<<std::endl;
      std::complex<FPType> res = std::tgamma(1-s)*std::pow(-w, s-1);
      constexpr FPType tp = 2.0 * M_PI;
      const std::complex<FPType> wup = w/tp;
      
      std::complex<FPType> w2 = wup;
      std::complex<FPType> pref = std::pow(tp, s)/M_PI;
      std::complex<FPType> gam = std::tgamma(1.0-s); //here we factor up the ratio of Gamma(1 - s + k)/k! . This ratio should be well behaved even for large k
      
      FPType sp, cp;
      sincos(M_PI/2.0 * s, &sp, &cp);
      /*Here we add the expression that would result from ignoring the zeta function in the series
       */
      std::complex<FPType> expis(cp, sp);
      std::complex<FPType> p = tp - std::complex<FPType>(0.0, 1.0) * w;
      std::complex<FPType> q = tp + std::complex<FPType>(0.0, 1.0) * w;
      res += std::complex<FPType>(0.0, 1.0) * gam * (conj(expis) * std::pow(p, s-1.0) - expis *std::pow(q, s-1.0));//this can be optimized for real values
      /*the above expression is the result of sum_k Gamma(1+k-s) /k! * sin(pi /2* (s-k)) * (w/2/pi)^k */
      /*therefore we only need to sample values of the zeta on the real axis that really differ from one*/
      res += pref * sp * gam * (std::tr1::__detail::__riemann_zeta(1-s) - 1.0);
      const unsigned int maxit = 200;
      unsigned int j = 1;
      bool terminate = false;
      gam*= (1.0 - s);
      while (!terminate)//assume uniform convergence
      {
	FPType rzarg = 1 + j - s;
	FPType rz = (std::tr1::__detail::__riemann_zeta(rzarg) - 1.0);//only the difference to one is needed
	FPType sine;
	if(j & 1)//save the reperated calculation of the sines
	{/*odd*/
	  sine = cp;
	  if ( !((j-1)/ 2 & 1) )
	    sine = -sine;
	}
	else
	{/*even*/
	  sine = sp;
	  if((j/2) & 1)
	    sine = -sine;
	}
	std::complex<FPType> nextterm =  w2 * gam * (sine * rz);
//	std::sin(M_PI/2.0*(s-j));
//	std::cout<<j<<" "<<nextterm<<" "<<rz<<std::endl;
	w2 *= wup;
	++j;
	gam  *= rzarg/(j);//equal to 1/(j+1) since we have incremented j in the line above
	terminate = (fpequal( std::abs(res + pref*nextterm), std::abs(res) ) || (j > maxit)) && !fpequal(std::abs(rz), 0.0);
	res += pref*nextterm;
      }
      std::cout<<"Iterations in PolyLogExp_neg: "<<j<<std::endl;
      return res;
    }
    static inline std::complex<FPType> PolyLog_Exp(const FPType s, std::complex<FPType> w)
    {
      /*reduce the imaginary part to the range where the series converge quickly*/
      if (fpequal<FPType>(std::rint(s), s))
      {//capture the cases of positive integer index
	int nu = static_cast<int> (lrint(s));
	std::cout<<nu<<std::endl;
	if(0 == nu)
	{
	  std::complex<FPType> t = std::exp(w);
	  return t/(1.0 - t);
	}
	else if (1 == nu)
	  return -std::log(1.0 - std::exp(w));
	else if (nu > 1)
	{//FIXME: check for real or non-real argument. asymptotic expansions
	  while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
	  while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
	    return std::tr1::__detail::__riemann_zeta(s);
	  else
	    return PolyLog_Exp(static_cast<uint>(nu) , w);
	}
	else
	  return PolyLog_Exp_neg(s, w);
      }
      else
      {//FIXME: check for real or non-real argument. asymptotic expansions
	  while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
	  while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	if (s < 0)
	  return PolyLog_Exp_neg(s, w);
	else
	{
	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
	    return std::tr1::__detail::__riemann_zeta(s);
	  else
	    return PolyLog_Exp_pos(s, w);
	}
      }
    }
    static inline std::complex<FPType> hurwitz_zeta(const FPType s, std::complex<FPType> x)
    {
      if ( 
	((x.imag() >= 0) && ((x.real() >= 0.0) && (x.real() <  1.0) )) || 
	((x.imag() <  0) && ((x.real() >  0.0) && (x.real() <= 1.0) ))
      )
      {
	constexpr FPType tp = 2.0*M_PI;
	FPType t = 1.0-s;
	std::complex<FPType> lpe = PolyLog_Exp(t, std::complex<FPType>(0.0, tp)/*== 2 pi I */ * x );
	//FIXME: This prefactor is prone to overflow
	return std::tgamma(t)* std::pow(tp, -t)* ( std::exp(std::complex<FPType>(0.0, -M_PI/2.0 * t)) * lpe + std::exp(std::complex<FPType>(0.0, M_PI/2.0 * t)) * conj(lpe) );
      }
      else 
      {
	std::cout<<"domain not supported!!"<<std::endl;
      }
    }
    static std::complex<FPType>* pl_store;
};

template <typename FPType>
FPType PowerLawRashbaChain<FPType>::disp(int k, FPType t, FPType my, FPType lambda) throw()//some system function seems to be also called epsilon
{
  FPType p = std::atan(lambda);
  return std::sqrt(t*t + lambda * lambda)* (exp(std::complex<FPType>(0.0,p)) * pl_store[k] ).real();
}

// template <typename FPType>
// FPType PowerLawRashbaChain<FPType>::GreenNullKSpacediag(FPType k, FPType tau, FPType t, FPType my, FPType lambda) throw()
// {
//   FPType ep = disp(k, t, my, lambda);
//   FPType em = disp(k, t, my, -lambda);
//   return 0.5*(std::exp( ep * tau) * fermi(beta * ep) + std::exp(em * tau) * fermi(beta * em));
// }
// 
// template <typename FPType>
// FPType PowerLawRashbaChain<FPType>::GreenNullKSpaceoffdiag(FPType k, FPType tau, FPType t, FPType my, FPType lambda) throw()
// {//remember to add the correct i or -i
//   FPType ep = disp(k, t, my, lambda);
//   FPType em = disp(k, t, my, -lambda);
//   return 0.5*(std::exp( ep * tau) * fermi(beta * ep) - std::exp(em * tau) * fermi(beta * em));
// }

template <typename FPType>
FPType PowerLawRashbaChain<FPType>::GreenNullKSpace(int k, FPType tau, FPType t, FPType my, FPType lambda) throw()
{
  FPType ep = disp(k, t, my, lambda);
  FPType retval;
//  std::cout<<ep*tau <<"> "<<std::log(0.01*std::numeric_limits<FPType>::max())<<std::endl;
  if(ep*tau > std::log(0.01*std::numeric_limits<FPType>::max()))
    retval = 0.5*std::exp(ep*(tau - beta/2))/std::cosh(ep*beta/2);
  else
    retval = std::exp( ep * tau) * fermi(beta * ep);
  return retval;
}

template<typename FPType>
std::complex<FPType>* PowerLawRashbaChain<FPType>::g = NULL;

template<typename FPType>
std::complex<FPType>* PowerLawRashbaChain<FPType>::pl_store = NULL;

template<typename FPType>
FPType PowerLawRashbaChain<FPType>::betaslice;

template<typename FPType>
FPType PowerLawRashbaChain<FPType>::beta;

template<typename FPType>
unsigned int PowerLawRashbaChain<FPType>::N;

template<typename FPType>
unsigned int PowerLawRashbaChain<FPType>::slices;

template <typename FPType>
void PowerLawRashbaChain<FPType>::tidyup()
{
    delete [] g;
}

template <typename FPType>
template <class CFG>
void PowerLawRashbaChain<FPType>::init(CFG& curparams)
{
    //lattice sites
    N = curparams.N;
    beta = curparams.beta;
    slices = 100000;//Number of TimeSlices
    const unsigned int slicesp = slices + 1;
    betaslice = beta / static_cast<FPType>(slices);
    pl_store = new std::complex<FPType>[curparams.N];
    for(uint k = 0; k < curparams.N; ++k)
    {
      pl_store[k] = PolyLog_Exp(2.0 + curparams.alpha,  std::complex<FPType>(0.0, (2*k)*M_PI/curparams.N) );
      std::cout<<k<<" -> "<<pl_store[k]<<std::endl;
    }
    g = new std::complex<FPType>[slicesp * N];
//     if (N == 1)//take special care of the 1 site hubbard model
//     {
//         for (unsigned int j = 0; j < slicesp; ++j)
//         {
//             g[j] = std::complex<FPType>(GreenNullKSpacediag(0, j * betaslice, curparams.t, curparams.mu, curparams.lambda),
//                                         GreenNullKSpaceoffdiag(0, j * betaslice, curparams.t, curparams.mu, curparams.lambda));
// //            std::cout<<j * betaslice<<" "<<g[j][0]<<std::endl;
//         }
//         return;
//     }
#pragma omp parallel for
    for (int j = 0; j < static_cast<int>(N); ++j)//for every realspacepoint
    {
      std::cout<<"j: "<<j<<std::endl;
        for (int i = 0; i < static_cast<int>(slicesp); ++i)//for every timeslice
        {
            std::complex<FPType> tempsum = 0;
            for (uint k = 0; k < N; ++k)//sum over all k-space values
            {
	      const std::complex<FPType> expkj = exp(std::complex<FPType>(0.0, 2.0 * M_PI/N * static_cast<FPType>(k * j)) );
	      double temp = GreenNullKSpace(k, i * betaslice, curparams.t, curparams.mu, curparams.lambda);
//	      std::cout<<k<<" "<<temp<<std::endl;
//	      if ( std::isnan(temp) ) exit(-1);
              tempsum +=  expkj * temp;//R_j = j, because the lattice constant is assumed to be equal to one.
	    }
            g[j*slicesp + i] = tempsum/static_cast<FPType>(N); 
        }
    }
    //Debugging Output
/*    std::ofstream dbg("dbg.txt");
    for(uint n = 0; n < N; ++n)
    {
        for(unsigned int i = 0; i < slicesp; ++i)
        {
          dbg<<(i * betaslice)<<" "<<-imag(g[n*slicesp + i]);
          dbg<<std::endl;
        }
        dbg<<"&"<<std::endl;
    }*/
//    exit(-1);
    delete [] pl_store;
    return;
}

template<typename FPType>
typename PowerLawRashbaChain<FPType>::FreeGreensFunctionReturnValueType PowerLawRashbaChain<FPType>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    //determine the Differences between the two
    FPType delta_tau = v1.tau - v2.tau;
    int delta = v1.site - v2.site;
    //Take care of negative values
    if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
    const std::complex<FPType> *const addr = g + delta*(slices + 1);//addr points to the memory area where the GF of the particular site resides.
    //Take care of negative values
    uint signchanges = 0;
    while (delta_tau < 0)
    {
        delta_tau += beta;
        ++signchanges;
    }
    while (delta_tau > beta)
    {
        delta_tau -= beta;
        ++signchanges;
    }
    FPType sign = (signchanges & 1? -1.0 : 1.0);
    if (fpequal(v1.tau, v2.tau))
    {
        //return only the particle number
        return sign * addr[0];
    }
    if (fpequal(delta_tau, beta))
        return sign * addr[slices];
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    long int tau_idx0 = lround(fptau_idx0);
    return lerp(rem, addr[tau_idx0], addr[tau_idx0 + 1]) * sign;//return the value of the greensfunction
}
#endif