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
#ifndef EXPONENTIAL_RASHBA_CHAIN_H
#define EXPONENTIAL_RASHBA_CHAIN_H
#include <complex>
#include <limits>
#include "ddqmc.h"
#include "Vertex.h"
#include "Greensfunction.h"
/**
An implementation of a 1D chain with an added Rashba-type interaction.
the transition rates decay exponentially
I assume that the lattice constant is a = 1
*/
template<typename FPType = float>
class ExponentialRashbaChain
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
    static inline FPType disp(FPType k, FPType t, FPType my, FPType lambda, FPType alpha) throw();
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
    static inline FPType GreenNullKSpace(FPType k, FPType tau, FPType t, FPType my, FPType lambda, FPType alpha) throw();
};

template <typename FPType>
FPType ExponentialRashbaChain<FPType>::disp(FPType k, FPType t, FPType my, FPType lambda, FPType alpha) throw()//some system function seems to be also called epsilon
{
  FPType p = std::atan(lambda);
  return std::sqrt(t*t + lambda * lambda)* std::cos(k + p)/(2*std::cosh(alpha) - 2.0*std::cos(k));
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
FPType ExponentialRashbaChain<FPType>::GreenNullKSpace(FPType k, FPType tau, FPType t, FPType my, FPType lambda, FPType alpha) throw()
{
  FPType ep = disp(k, t, my, lambda, alpha);
  FPType retval;
//  std::cout<<ep*tau <<"> "<<std::log(0.01*std::numeric_limits<FPType>::max())<<std::endl;
  if(ep*tau > std::log(0.01*std::numeric_limits<FPType>::max()))
    retval = 0.5*std::exp(ep*(tau - beta/2))/std::cosh(ep*beta/2);
  else
    retval = std::exp( ep * tau) * fermi(beta * ep);
  return retval;
}

template<typename FPType>
std::complex<FPType>* ExponentialRashbaChain<FPType>::g = NULL;

template<typename FPType>
FPType ExponentialRashbaChain<FPType>::betaslice;

template<typename FPType>
FPType ExponentialRashbaChain<FPType>::beta;

template<typename FPType>
unsigned int ExponentialRashbaChain<FPType>::N;

template<typename FPType>
unsigned int ExponentialRashbaChain<FPType>::slices;

template <typename FPType>
void ExponentialRashbaChain<FPType>::tidyup()
{
    delete [] g;
}

template <typename FPType>
template <class CFG>
void ExponentialRashbaChain<FPType>::init(CFG& curparams)
{
    //lattice sites
    N = curparams.N;
    beta = curparams.beta;
    slices = 100000;//Number of TimeSlices
    const unsigned int slicesp = slices + 1;
    betaslice = beta / static_cast<FPType>(slices);
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
        for (int i = 0; i < static_cast<int>(slicesp); ++i)//for every timeslice
        {
            std::complex<FPType> tempsum = 0;
            for (uint k = 0; k < N; ++k)//sum over all k-space values
            {
	      const std::complex<FPType> expkj = exp(std::complex<FPType>(0.0, 2.0 * M_PI/N * static_cast<FPType>(k * j)) );
	      FPType temp = GreenNullKSpace(2*k*M_PI/N, i * betaslice, curparams.t, curparams.mu, curparams.lambda, curparams.alpha);
//	      std::cout<<k<<" "<<temp<<std::endl;
//	      if ( std::isnan(temp) ) exit(-1);
              tempsum +=  expkj * temp;//R_j = j, because the lattice constant is assumed to be equal to one.
	    }
            g[j*slicesp + i] = tempsum/static_cast<FPType>(N); 
        }
    }
    //Debugging Output
    /*std::ofstream dbg("dbg.txt");
    for(uint n = 0; n < N; ++n)
    {
        for(unsigned int i = 0; i < slicesp; ++i)
        {
          dbg<<(i * betaslice)<<" "<<-imag(g[n*slicesp + i]);
          dbg<<std::endl;
        }
        dbg<<"&"<<std::endl;
    }
    exit(-1);*/
    return;
}

template<typename FPType>
typename ExponentialRashbaChain<FPType>::FreeGreensFunctionReturnValueType ExponentialRashbaChain<FPType>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    //determine the differences between the two
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