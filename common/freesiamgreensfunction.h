/*
    Copyright (c) 2010, Florian Goth
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

#ifndef FREESIAMGREENSFUNCTION_H
#define FREESIAMGREENSFUNCTION_H

#include <complex>
#include "Greensfunction.h"
#include "Vertex.h"
#include <limits>
/**
This class is the Greensfunction for the Single Impurity Anderson model(SIAM) in the Wide-band Limit
$W->\infty$ under the condition V/W = const.
*/

template<typename FPType_ = float>
class G_SIAM_Omega
{
  public:
    typedef FPType_ FPType;
    typedef std::complex<FPType> RetType;
    enum {
      has_real_FourierTransform = true,
      has_TRS = true,
    };
    G_SIAM_Omega(FPType v_, FPType ed_, FPType w_, FPType my_) : v(v_), ed(ed_), w(w_), my(my_) {}
    RetType operator()(FPType omegan)
    {
 //FPType phi = atan(2.0*omegam / w);
/*
    return static_cast<FPType>(1.0)/(std::complex<FPType>(-ed, -omegam ) + 
    v*v/w *log( (std::complex<FPType>(0.0, 1.0) * omegam + 0.5*w)/(std::complex<FPType>(0.0, 1.0) * omegam - 0.5*w))  );
  */
  /*  
if (omegam > 0)//the branching is necessary to restrict the value of PI+2*phi to the range of the arg() which is [-PI, PI]
    return static_cast<FPType>(1.0)/std::complex<FPType>(-ed, -omegam + v*v/w*(-M_PI + 2.0*phi));
    else
    return static_cast<FPType>(1.0)/std::complex<FPType>(-ed, -omegam + v*v/w*(+M_PI + 2.0*phi));
*/
if(omegan > 0)//Wide Band Limit
      return static_cast<FPType>(1.0)/std::complex<FPType>(-ed, -omegan - v*v/w*M_PI);
      else
      return static_cast<FPType>(1.0)/std::complex<FPType>(-ed, -omegan + v*v/w*M_PI);
    }
  private:
    FPType v;
    FPType ed;
    FPType w;
    FPType my;
};

template<typename FPType = float>
class FreeSIAMGreensFunction
{
public:
    enum {timeevolution = 0,
    has_Spin_symmetry = true,
    has_Giomega = true,
    Is_Impurity_model = true,
    has_TRS = true,
    };
    typedef Basic_Vertex<FPType> Vertex;
    typedef FPType FreeGreensFunctionReturnValueType;///< a typedef of the data type of the free greensfunction
    typedef G_SIAM_Omega<FPType> GOmega;
    typedef typename GOmega::RetType GOmegaRetType;
    /**
    This evaluates the value of the free particle Greens-function for the two given vertices
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
    To access the nr of atoms.(In the SIAM this is mostly for compatability)
    @return the nr of atoms in the chain
    */
    static inline unsigned int getLen() throw()
    {
        return 1;
    }
    /**
    Get the length of the interactioninterval.
    @return the inverse temperature beta
    */
    static FPType getContourLen() throw()
    {
        return beta;
    }
    /**
    The free particle Greens function in Matsubara frequencies
    */
    inline static typename G_SIAM_Omega<FPType>::RetType gomega(FPType om, Vertex& v1, Vertex& v2)
    {
      if(v1.spin == v2.spin)
      return (*gomega_int)(om);
      else
	return 0.0;
    }
private:
    static G_SIAM_Omega<FPType >* gomega_int;
    static FreeGreensFunctionReturnValueType* g;///< the free particle greensfunction
    static unsigned int slices;///< the number of timeslices. this means the resolution of the free greensfunction in tau space
    static FPType betaslice;///<the length of one timeslice on the tau axis
    static FPType beta;///<the inverse temperature beta
    static FPType v;
    static FPType w;
    static FPType ed;
    static FPType my;
};

template<typename FPType>
typename FreeSIAMGreensFunction<FPType>::FreeGreensFunctionReturnValueType* FreeSIAMGreensFunction<FPType>::g = NULL;

template <typename FPType>
G_SIAM_Omega<FPType>* FreeSIAMGreensFunction<FPType>::gomega_int = NULL;

template<typename FPType>
FPType FreeSIAMGreensFunction<FPType>::betaslice;

template<typename FPType>
FPType FreeSIAMGreensFunction<FPType>::beta;

template<typename FPType>
unsigned int FreeSIAMGreensFunction<FPType>::slices;

template <typename FPType>
void FreeSIAMGreensFunction<FPType>::tidyup()
{
    delete [] g;
    delete gomega_int;
}

template <typename FPType>
template <class CFG>
void FreeSIAMGreensFunction<FPType>::init(CFG& curparams)
{
    //lattice sites
    beta = curparams.beta;
    slices = 10000;//Number of TimeSlices
    const unsigned int slicesp = slices + 1;
    betaslice = beta / static_cast<FPType>(slices);
    gomega_int = new G_SIAM_Omega<FPType>(curparams.V, curparams.ed, curparams.W, curparams.mu);
    g = new FreeGreensFunctionReturnValueType[slicesp];
    matsubarafouriertransform(*gomega_int, g, beta, betaslice, slicesp, slices);
    return;
}

template<typename FPType>
typename FreeSIAMGreensFunction<FPType>::FreeGreensFunctionReturnValueType FreeSIAMGreensFunction<FPType>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    const FPType tiny = std::numeric_limits<FPType>::epsilon();
    //determine the Differences between the two
    FPType delta_tau = v1.tau - v2.tau;
    //Take care of negative values
    if (std::abs(delta_tau) < tiny)
    {
        //return only the particle number
        return g[0];
    }
    if (std::abs(delta_tau - beta) < tiny)
        return -g[slices];
    FPType sign = 1.0;
    if (delta_tau < 0)//Again take care of negative values
    {
        delta_tau += beta;
        sign = -1.0;//The Matsubara GreensFunction changes sign if we add a beta!!!!!!!!!!
    }
    if (delta_tau < 0)//Take care of very negative values
    {
        delta_tau += beta;
        sign *= -1.0;//The Matsubara GreensFunction changes sign if we add a beta!!!!!!!!!!
    }
    if(delta_tau > beta)
    {
      delta_tau -= beta;
      sign *= -1.0;
    }
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    long int tau_idx0 = lround(fptau_idx0);
//    std::cout<<"tau_0: "<<tau_idx0<<" "<<g[tau_idx0]<<std::endl;
    return lerp(rem, g[tau_idx0], g[tau_idx0 + 1]) * sign;//return the value of the greensfunction
}

#endif // FREESIAMGREENSFUNCTION_H
