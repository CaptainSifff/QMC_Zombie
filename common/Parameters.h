/***************************************************************************
 *   Copyright (C) 2009-2013 by Florian Goth   *
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
#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <inttypes.h>
#include <cmath>
#include <iostream>
#include <valarray>
#include "Models.h"
#define DP(text) std::cout << (#text) << " = " << text << std::endl;

struct Parameters
{
    Models model;
    uint32_t N;
    double t;
    double U;
    double beta;
    double mu;
    double t_exp;
    double delta;
    double V;
    double W;
    double ed;
    
    uint32_t Nx;
    uint32_t Nb;
    double lambda;
    double alpha;
    uint32_t N_wait;
    uint32_t prebins;///< Here we store the number of measurements over which to average before sending to the server.
    double contourlen;
    double delta_s;
    unsigned int functionpoints;
    char* data;
    uint datalen;
    template <class Net>
    Parameters(Net& z, Models model) : data(NULL)
    {
      switch(model)
      {
	case HUBBARD_CHAIN:
	case COLD_ATOMS:
	    N = z.template recvblock<uint32_t>();
            DP(N);
            V = -100000;
            W = -100000;
            ed = -100000;
	  break;
	case SIAM:
	    //FIXME probably introduce booleans like dot-model and lattice-model for proper propagation of the lattice length
            N = 1;
            V = z.template recvblock<double>();
            DP(V);
            W = z.template recvblock<double>();
            DP(W);
            ed = z.template recvblock<double>();
            DP(ed);
	  break;
	case IMAG_MODEL_FROM_FILE:
	{
	    N = z.template recvblock<uint32_t>();
            DP(N);
            V = -100000;
            W = -100000;
            ed = -100000;
	    std::valarray<char> temp(z.template recvblock<std::valarray<char> >());
	    data = new char[temp.size()];
	    for(uint k = 0; k < temp.size(); ++k)
	      data[k] = temp[k];
	    datalen = static_cast<uint>(temp.size());
	  break;
	}
	case KONDO_IMP_TI:
	    //FIXME probably introduce booleans like dot-model and lattice-model for proper propagation of the lattice length
            N = 1;
            V = z.template recvblock<double>();
            DP(V);
            ed = z.template recvblock<double>();
	    DP(ed);
	    Nx = z.template recvblock<uint32_t>();
	    DP(Nx);
	    Nb = z.template recvblock<uint32_t>();
	    DP(Nb);
	    lambda = z.template recvblock<double>();
	    DP(lambda);
	  break;
	case RASHBA_CHAIN:
	    N = z.template recvblock<uint32_t>();
            DP(N);
            V = -100000;
            W = -100000;
            ed = -100000;
	    lambda = z.template recvblock<double>();
	    DP(lambda);
	  break;
	case RASHBA_CHAIN_EXPONENTIAL:
	case RASHBA_CHAIN_POWER_LAW:
	    N = z.template recvblock<uint32_t>();
            DP(N);
            V = -100000;
            W = -100000;
            ed = -100000;
	    lambda = z.template recvblock<double>();
	    alpha = z.template recvblock<double>();
	    DP(lambda);
	    DP(alpha);
	  break;
      }
        t = z.template recvblock<double>();
        DP(t);
        U = z.template recvblock<double>();
        DP(U);
        beta = z.template recvblock<double>();
        DP(beta);
        mu = z.template recvblock<double>();
	std::cout.precision(10);
	std::cout<<mu<<std::endl;
        DP(mu);
        t_exp = z.template recvblock<double>();
        DP(t_exp);
        delta = z.template recvblock<double>();
        DP(delta);
        N_wait = z.template recvblock<uint32_t>();
	prebins = z.template recvblock<uint32_t>();
        DP(N_wait);
	DP(prebins);
        delta_s = z.template recvblock<double>();
        DP(delta_s);
        contourlen = 2.0 * t_exp + beta;
        functionpoints = static_cast<unsigned int>(std::round(t_exp == 0. ? beta/delta_s : t_exp/delta_s));
    }
};
#endif
