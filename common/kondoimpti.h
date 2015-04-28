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

#ifndef KONDOIMPTI_H
#define KONDOIMPTI_H

#include <complex>
#include <fstream>
#include "Greensfunction.h"
#include "Vertex.h"
#include "MTL/Matrix.h"
#include <limits>
#include "MTL/CommonTypes.h"
/**
some class for the Kondo Impurity in a 2D TI bath
*/

template <typename FPType>
struct GOmegaData
{
    FPType lambda;
    std::complex<FPType> u;
    std::complex<FPType>* evec;
};

template<typename FPType_ = float>
class G_TI_Omega
{
public:
    typedef FPType_ FPType;
    typedef std::complex<FPType> RetType;
    typedef RetType value_type;
    enum {
        has_real_FourierTransform = false
    };
    G_TI_Omega(const GOmegaData<FPType> *const mydat, uint l, uint ks, FPType beta, FPType v_, FPType ed_) : data(mydat), v(v_), ed(ed_), betaoverpi(beta / M_PI), ksize(ks), len(l)
    {
        gomegaplus = new RetType[omega_max];
        gomegaminus = new RetType[omega_max];
#pragma omp parallel for
        for (uint m = 0; m < omega_max; ++m)
        {
            std::complex<FPType> omegam = std::complex<FPType>(0.0, (2*m+1)*M_PI/beta);
            RetType sumplus = 0.0;
            RetType summinus = 0.0;
            for (uint k = 0; k < len; ++k)
                for (uint n = 0; n < ksize; ++n)
                {
                    sumplus += norm(data[k*ksize + n].u)/(omegam - data[k*ksize + n].lambda);
                    summinus += norm(data[k*ksize + n].u)/(-omegam - data[k*ksize + n].lambda);
                }
            sumplus /= len;
            summinus /= len;
            gomegaplus[m] = 1.0/(-omegam - ed - v*v*conj(sumplus));
            gomegaminus[m] = 1.0/(omegam - ed - v*v*conj(summinus));
        }
        return;
    }
    inline RetType operator()(FPType omegan) const
    {
        RetType* g = gomegaplus;
        if (omegan < 0)
        {
            omegan = -omegan;
            g = gomegaminus;
        }
        long int idx = lround(omegan*betaoverpi);
        --idx;
        idx/=2;
        if (idx < omega_max)
            return g[idx];
        else return RetType (0.0, 1.0/omegan);
        /*      if(omegan > 0)//Wide Band Limit
              return static_cast<FPType>(1.0)/std::complex<FPType>(-ed, -omegan - 0.75*0.75/4.0*M_PI);
              else
              return static_cast<FPType>(1.0)/std::complex<FPType>(-ed, -omegan + 0.75*0.75/4.0*M_PI);*/
    }
    inline RetType operator()(int omegan) const
    {
        RetType* g = gomegaplus;
        if (omegan < 0)
        {
            omegan = -omegan-1;
            g = gomegaminus;
        }
        if (omegan < static_cast<int>(omega_max))
            return g[omegan];
        else return RetType (0.0, 1.0/((2*omegan + 1)/betaoverpi));
    }
    ~G_TI_Omega()
    {
        delete [] gomegaplus;
        delete [] gomegaminus;
	for (uint k = 0; k < len; ++k)
                for (uint n = 0; n < ksize; ++n)
		  delete [] data[k*ksize + n].evec;
	delete [] data;
    }
    uint kpoints() const {return ksize;}
    const GOmegaData<FPType> *const data;
private:
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (GCC_VERSION > GCC_VER(4,5,0))
    static constexpr uint omega_max = 40000;
#else
    static const uint omega_max = 40000;
#endif
    FPType v;
    FPType ed;
    RetType* gomegaplus;
    RetType* gomegaminus;
    FPType betaoverpi;
    uint ksize;
    uint len;
};

template<typename FPType = float>
class KondoImpTI
{
public:
    enum {timeevolution = 0,
          has_Spin_symmetry = false,
          has_Giomega = true,
          Is_Impurity_model = true,
	  has_TRS = true,
         };
    typedef Basic_Vertex<FPType> Vertex;
    typedef std::complex<FPType> FreeGreensFunctionReturnValueType;///< a typedef of the data type of the free greensfunction
    typedef G_TI_Omega<FPType> GOmega;
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
    inline static typename G_SIAM_Omega<FPType>::RetType gomega(FPType om, Vertex& v1, Vertex& v2) throw()
    {
      if (v1.spin != v2.spin) return 0.0;
      else
        return v1.spin == UP ? (*gomegaup)(om) : (*gomegadown)(om);
    }
    static G_TI_Omega<FPType>* gomegaup;
    static G_TI_Omega<FPType>* gomegadown;
    static FPType gethybridization() {return v;}
private:
    static FreeGreensFunctionReturnValueType** g;///< the free particle greensfunction
    static unsigned int slices;///< the number of timeslices. this means the resolution of the free greensfunction in tau space
    static FPType betaslice;///<the length of one timeslice on the tau axis
    static FPType beta;///<the inverse temperature beta
    static FPType ed;///< dot-level
    static FPType my;///< chemical potential
    static FPType v;///< hybridization
};

template<typename FPType>
typename KondoImpTI<FPType>::FreeGreensFunctionReturnValueType** KondoImpTI<FPType>::g = NULL;

template <typename FPType>
G_TI_Omega<FPType>* KondoImpTI<FPType>::gomegaup = NULL;

template <typename FPType>
G_TI_Omega<FPType>* KondoImpTI<FPType>::gomegadown = NULL;

template<typename FPType>
FPType KondoImpTI<FPType>::betaslice;

template<typename FPType>
FPType KondoImpTI<FPType>::v;

template<typename FPType>
FPType KondoImpTI<FPType>::ed;

template<typename FPType>
FPType KondoImpTI<FPType>::beta;

template<typename FPType>
unsigned int KondoImpTI<FPType>::slices;

template <typename FPType>
void KondoImpTI<FPType>::tidyup()
{
    delete [] g[0];
    delete [] g[1];
    delete [] g;
    delete gomegaup;
    delete gomegadown;
}

template <typename FPType>
template <class CFG>
void KondoImpTI<FPType>::init(CFG& curparams)
{
    v = curparams.V;
    ed = curparams.ed;
    uint Nx = curparams.Nx;
    uint Nb = curparams.Nb;
    uint bsize = 2;
    uint ksize = Nb * bsize;
    uint size = Nx * ksize;
    FPType lambda = curparams.lambda;
    beta = curparams.beta;
    slices = 10000;//Number of TimeSlices
    const unsigned int slicesp = slices + 1;
    betaslice = beta / static_cast<FPType>(slices);
    g = new FreeGreensFunctionReturnValueType*[2];
    G_TI_Omega<FPType>* gomegaarr[2];
    std::ofstream summen("summen.txt");
    std::ofstream gdownfile("gdown.txt");
    std::ofstream gupfile("gup.txt");
    std::ofstream energies("energies.txt");
    for (int s = 0; s < 2; ++s)
    {
        GOmegaData<FPType>* data = new GOmegaData<FPType>[Nx*ksize];
        g[s] = new FreeGreensFunctionReturnValueType[slicesp];
        int sigma = 2*s - 1;
        for (uint kidx = 0; kidx < Nx; ++kidx)
        {
            FPType k = kidx*2.0*M_PI/Nx;
            summen<<k<<std::endl;
            std::complex<FPType> expkhalf = std::exp(std::complex<FPType>(0.0, kidx*M_PI/Nx));
	    std::complex<FPType> expk = exp(std::complex<FPType>(0.0, k));
	    std::complex<FPType> cexpk = conj(expk);
	    FPType ls = -2.0 *lambda*sigma;
	    FPType sinklambda = imag(expk)* ls;
	    FPType cosfack = 2.0*std::real(expkhalf)/*cos(kidx*M_PI/Nx)*/;
	    FPType sinfack = std::imag(expkhalf);
	    std::complex<FPType> lambdaentry = - sinfack * conj(expkhalf)*ls;
	    std::complex<FPType> clambdaentry = conj(lambdaentry);
	    std::complex<FPType> entry = cosfack*expkhalf;//= 1 + expk
	    std::complex<FPType> centry = conj(entry);//cosfack*cexpk;// = 1 + cexpk
            CmplxMat hzero_k_sigma(ksize, ksize, MTLICCL::Identity<FPType>(0.0, 0.0));
//            CmplxMat hlambda_k_sigma(ksize, ksize, MTLICCL::Identity<FPType>(0.0, 0.0));
            for (uint b = 0; b < Nb; ++b)
            {
//	      hzero_k_sigma(b*bsize + 0, b *bsize + 0) = 0.1;
//	      hzero_k_sigma(b*bsize + 1, b *bsize + 1) = 0.1;

                 hzero_k_sigma(b*bsize + 1, b *bsize + 0) = -centry;
                 hzero_k_sigma(b*bsize + 0, b *bsize + 1) = -entry;
                if (b != 0)
                {
                    hzero_k_sigma((b-1)*bsize + 1, b *bsize + 0) = -cexpk;
                    hzero_k_sigma(b*bsize + 0, (b-1) *bsize + 1) = -expk;
		    
		    //Now the Spin Orbit parts
		    hzero_k_sigma((b-1)*bsize + 0, b *bsize + 0) = lambdaentry;
                    hzero_k_sigma(b*bsize + 0, (b-1) *bsize + 0) = clambdaentry;
                    
                    hzero_k_sigma((b-1)*bsize + 1, b *bsize + 1) = -lambdaentry;
                    hzero_k_sigma(b*bsize + 1, (b-1) *bsize + 1) = -clambdaentry;
                }
                 hzero_k_sigma(b * bsize + 0, b * bsize + 0) = sinklambda;
                 hzero_k_sigma(b * bsize + 1, b * bsize + 1) = -sinklambda;
            }
//            std::cout<<(std::complex<double>(lambda*sigma,0.0)*hlambda_k_sigma-hzero_k_sigma)<<std::endl;
//            CmplxMat fullmat(/*std::complex<double>(lambda*sigma,0.0)*hlambda_k_sigma -*/ hzero_k_sigma);
            MTLICCL::EigenSystemSolver<CmplxMat> ess2(hzero_k_sigma);
            ess2.calculateEVs();
            CmplxMat zcplx(ess2.tridiagonalize());
            std::vector<MTLICCL::SolutionPair<double> > esystemcplx = ess2.eigensystem();
            for (uint n = 0; n < ksize; ++n)
            {
                data[kidx * ksize + n].lambda = esystemcplx[n].evalue;
		data[kidx * ksize + n].evec = new std::complex<FPType>[ksize];
                energies<<k<<" "<<esystemcplx[n].evalue<<std::endl;
                data[kidx * ksize + n].u = 0.0;
                for (uint q = 0; q < ksize; ++q)
		{
                    data[kidx * ksize + n].u += zcplx(0, q)*esystemcplx[n].evector[q];
		}
		for(uint i = 0; i < ksize; ++i)
		  for(uint q = 0; q < ksize; ++q)
		  {
		    data[kidx * ksize + n].evec[i] += zcplx(i, q)*esystemcplx[n].evector[q];
		  }
            }
            
            
            CmplxMat evscplx(zcplx.Rows(), zcplx.Rows());
            for (uint i = 0; i < zcplx.Rows(); ++i)
                for (uint j = 0; j < zcplx.Rows(); ++j)
                {
                    evscplx(j, i) = esystemcplx[i].evector[j];
                }
//   std::cout<<zcplx*evscplx<<std::endl;
            CmplxMat zhc(~(zcplx*evscplx));
            for (uint j = 0; j < zhc.Rows(); ++j)
                for (uint i = 0; i < zhc.Columns(); ++i)
                    zhc(j,i) = conj(zhc(j,i));
            CmplxMat id(zhc*(zcplx*evscplx));
            for (uint j = 0; j < id.Rows(); ++j)
            {
                double sum = 0.0;
                for (uint i = 0; i < id.Rows(); ++i) {
                    if (i!=j)
                        sum += abs(id(i,j));
                }
                std::cout.precision(12);
                summen<<"Zeilensumme "<<j<<"    : "<<std::scientific<<sum<<" + "<<std::abs(id(j,j))<<std::endl;
            }
        }
        energies<<"&"<<std::endl;
        std::cout<<"finished k-sum"<<std::endl;
        gomegaarr[s] = new G_TI_Omega<FPType>(data, Nx, ksize, beta, curparams.V, curparams.ed);
        std::cout<<"finished gomega"<<std::endl;
        matsubarafouriertransform(*(gomegaarr[s]), g[s], beta, betaslice, slicesp, slices);
        for (uint k = 0; k < slicesp; ++k)
            (s == 0? gdownfile: gupfile)<<k*betaslice<<" "<<std::scientific<<real(g[s][k])<<std::endl;
//	exit(-1);
    }
    gomegadown = gomegaarr[0];
    gomegaup = gomegaarr[1];
    return;
}

template<typename FPType>
typename KondoImpTI<FPType>::FreeGreensFunctionReturnValueType KondoImpTI<FPType>::eval(const Vertex& v1, const Vertex& v2) throw()
{
    const FPType tiny = std::numeric_limits<FPType>::epsilon();
    //determine the Differences between the two
    FPType delta_tau = v1.tau - v2.tau;
    if (v1.spin != v2.spin) return 0.0;
    const FreeGreensFunctionReturnValueType *const glocal = (v1.spin == UP ? g[1] : g[0]);
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
    if (std::abs(delta_tau) < tiny)
    {
        //return only the particle number
        return sign*glocal[0];
    }
    if (std::abs(delta_tau - beta) < tiny)
        return sign*glocal[slices];
    FPType fptau_idx0;
    FPType rem = std::modf(delta_tau/betaslice, &fptau_idx0);//round to the smaller index and determine how far we're of
    long int tau_idx0 = lround(fptau_idx0);
    return sign * glocal[tau_idx0];
//    return lerp(rem, glocal[tau_idx0], glocal[tau_idx0 + 1]) * sign;//return the value of the greensfunction
}
#endif // KONDOIMPTI_H
