/*
    Copyright (c) 2011, Florian Goth
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef IMAG_MODEL_FROM_SERVER_H
#define IMAG_MODEL_FROM_SERVER_H

#include <complex>
#include "Greensfunction.h"
#include "Vertex.h"
#include "charconverters.h"
#include <limits>
/**
A class for simulating Hubbard models with a Server supplied free Greensfunction
*/
template<typename T, bool Spin_Symmetry>
class ImagModelFromServer
{
public:
    enum {timeevolution = 0,
    has_Spin_symmetry = Spin_Symmetry,
    has_Giomega = false,
    Is_Impurity_model = false,
    has_TRS = false,
    };
    typedef Hubbard_Vertex<double> Vertex;
    typedef double FPType;//assume that we always get double data from the server
    typedef T FreeGreensFunctionReturnValueType;///< a typedef of the data type of the free greensfunction
    /**
    This evaluates the value of the free Greens-function for the two given vertices
    @param v1 the first vertex
    @param v2 the second vertex
    @return the value of the free greensfunction evaluated with the given vertices
    */
    static inline FreeGreensFunctionReturnValueType eval(const Vertex& v1, const Vertex& v2) throw();
    /**
    A function for initializing the tables that make up the Greensfunction. the necessary parameters are read from the parameter structs.
    */
    template <class CFG>
    static inline void init(CFG&);
    /**
    This function frees the memory used by the tables for the Greensfunction
    */
    static inline void tidyup();
    /**
    To access the nr of atoms.
    @return the nr of atoms in the chain
    */
    static inline unsigned int getLen() throw()
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
    static FreeGreensFunctionReturnValueType* g;///< the free particle greensfunction
    static unsigned int slices;///< the number of timeslices. this means the resolution of the free greensfunction in tau space
    static FPType betaslice;///<the length of one timeslice on the tau axis
    static FPType beta;///<the inverse temperature beta
    static FPType my;
    static unsigned int N;
    static unsigned int Nslices;
};

template<typename T, bool Spin_Symmetry>
typename ImagModelFromServer<T, Spin_Symmetry>::FreeGreensFunctionReturnValueType* ImagModelFromServer<T,Spin_Symmetry>::g = NULL;

template<typename T, bool Spin_Symmetry>
unsigned int ImagModelFromServer<T, Spin_Symmetry>::N;

template<typename T, bool Spin_Symmetry>
unsigned int ImagModelFromServer<T, Spin_Symmetry>::Nslices;

template<typename T, bool Spin_Symmetry>
double ImagModelFromServer<T, Spin_Symmetry>::betaslice;

template<typename T, bool Spin_Symmetry>
double ImagModelFromServer<T, Spin_Symmetry>::beta;

template<typename T, bool Spin_Symmetry>
unsigned int ImagModelFromServer<T, Spin_Symmetry>::slices;

template<typename T, bool Spin_Symmetry>
void ImagModelFromServer<T, Spin_Symmetry>::tidyup()
{
    delete [] g;
}

template<typename T, bool Spin_Symmetry>
template <class CFG>
void ImagModelFromServer<T, Spin_Symmetry>::init(CFG& curparams)
{
    //lattice sites
    beta = curparams.beta;
    N = curparams.N;
    g = new T[curparams.datalen/sizeof(T)];
    for(uint k = 0; k < curparams.datalen/sizeof(T); ++k)
    {
      g[k] = CH<T>(curparams.data + k*sizeof(T));
    }
    delete [] curparams.data;
    
    slices = curparams.datalen/sizeof(T)/N;
    if(!Spin_Symmetry) slices /= 2;
    Nslices = N*slices;
    betaslice = beta / static_cast<FPType>(slices);
    return;
}

/**
 * A word of caution: usually the spin in the vertex denotes the Ising Spin of a configuration
 * But just before we evaluate a specific greensfunction(hence, we call eval() )
 * The spin of a Vertex gets abused to denote the true fermion spin.(See create_All_four_Spin_Possiblilities()  ).
 */

//Now follow the full specializations that give the behaviour when evaluating the greensfunction
template<>
double ImagModelFromServer<double, true>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    const FPType tiny = betaslice;
    //determine the difference between the two
    FPType delta_tau = v1.tau - v2.tau;
    int delta = v1.site - v2.site;
    //Take care of negative values
    if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
    if (std::abs(delta_tau) < 0.5*tiny)
    {
        //return only the particle number
        return g[delta*slices];
    }
    if (delta_tau == beta)
        return g[(delta + 1) * slices - 1];
    FPType sign = 1.0;
    if (delta_tau < 0)//Again take care of negative values
    {
        delta_tau = beta + delta_tau;
        sign = -1.0;//The Matsubara GreensFunction changes sign if we add a beta!!!!!!!!!!
    }
    if(delta_tau > beta)
    {
        delta_tau -= beta;
        sign *= -1.0;
    }
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    std::size_t tau_idx0 = static_cast<std::size_t>(lround(fptau_idx0));
    return lerp(rem, g[delta*slices + tau_idx0], g[delta*slices + tau_idx0 + 1]) * sign;//return the value of the greensfunction
}

template<>
double ImagModelFromServer<double, false>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    if(v1.spin != v2.spin) return 0.0;
    uint baseindex = 0;
    if(v1.spin == DOWN) baseindex = Nslices;
    const FPType tiny = betaslice;
    //determine the difference between the two
    FPType delta_tau = v1.tau - v2.tau;
    int delta = v1.site - v2.site;
    //Take care of negative values
    if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
    if (std::abs(delta_tau) < 0.5*tiny)
    {
        //return only the particle number
        return g[baseindex + delta*slices];
    }
    if (delta_tau == beta)
        return g[baseindex + (delta + 1) * slices - 1];
    FPType sign = 1.0;
    if (delta_tau < 0)//Again take care of negative values
    {
        delta_tau = beta + delta_tau;
        sign = -1.0;//The Matsubara GreensFunction changes sign if we add a beta!!!!!!!!!!
    }
    if(delta_tau > beta)
    {
        delta_tau -= beta;
        sign *= -1.0;
    }
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    std::size_t tau_idx0 = static_cast<std::size_t>(lround(fptau_idx0));
    return lerp(rem, g[baseindex + delta*slices + tau_idx0], g[baseindex + delta*slices + tau_idx0 + 1]) * sign;//return the value of the greensfunction
}

template<>
std::complex<double> ImagModelFromServer<std::complex<double>, true>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    const FPType tiny = betaslice;
    //determine the difference between the two
    FPType delta_tau = v1.tau - v2.tau;
    int delta = v1.site - v2.site;
    //Take care of negative values
    if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
    if (std::abs(delta_tau) < 0.5*tiny)
    {
        //return only the particle number
        return g[delta*slices];
    }
    if (delta_tau == beta)
        return g[(delta + 1) * slices - 1];
    FPType sign = 1.0;
    if (delta_tau < 0)//Again take care of negative values
    {
        delta_tau = beta + delta_tau;
        sign = -1.0;//The Matsubara GreensFunction changes sign if we add a beta!!!!!!!!!!
    }
    if(delta_tau > beta)
    {
    delta_tau -= beta;
    sign *= -1.0;
    }
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    std::size_t tau_idx0 = static_cast<std::size_t>(lround(fptau_idx0));
    return g[delta*slices + tau_idx0] * sign;//return the value of the greensfunction
}

template<>
std::complex<double> ImagModelFromServer<std::complex<double>, false>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    if(v1.spin != v2.spin) return 0.0;
    uint baseindex = 0;
    if(v1.spin == DOWN) baseindex = Nslices;
    const FPType tiny = std::numeric_limits<FPType>::epsilon();
    //determine the difference between the two
    FPType delta_tau = v1.tau - v2.tau;
    int delta = v1.site - v2.site;
    //Take care of negative values
    if (delta < 0) delta = N + delta;//because delta is in this case already a negative number
    if (std::abs(delta_tau) < tiny)
    {
        //return only the particle number
        return g[baseindex + delta*slices];
    }
    if (std::abs(delta_tau - beta) < tiny)
        return -g[baseindex + (delta + 1) * slices - 1];
    FPType sign = 1.0;
    if (delta_tau < 0)//Again take care of negative values
    {
        delta_tau = beta + delta_tau;
        sign = -1.0;//The Matsubara GreensFunction changes sign if we add a beta!!!!!!!!!!
    }
    if(delta_tau > beta)
    {
    delta_tau -= beta;
    sign *= -1.0;
    }
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    std::size_t tau_idx0 = static_cast<std::size_t>(lround(fptau_idx0));
    return g[baseindex + delta*slices + tau_idx0] * sign;//return the value of the greensfunction
}
#endif
