/***************************************************************************
 *   Copyright (C) 2009 by Florian Goth   *
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//generally useful stuff
#include <iostream>
#include <cstdlib>

#ifdef __MINGW32__
#include <windows.h>
#define sleep(n) Sleep(1000 * n)
#endif
#include <unistd.h>
#include <ctime>
//stuff for measuring the amount of allocations in the matrices
unsigned long memctr = 0;
unsigned long memsize = 0;
static clock_t start = clock();
static clock_t endtime;
//stuff for network
#include "ParallelClient.h"
#include "commands.h"

//stuff for the simulation
#include "prng.h"
#include "sampling_algorithms.h"
#include "FreeEquilibriumsGreensFunction.h"
#include "freesiamgreensfunction.h"
#include "imag_model_from_server.h"
#include "kondoimpti.h"
#include "rashba_chain.h"
#include "powerlawRashba.h"
#include "exponentialRashba.h"

#include "QMC_Main.h"
#include <vector>
#include <stdexcept>
#include "Models.h"

using namespace std;

template <typename FPType, class DDConfSIAM>
struct DDQMC_Init_Helper<FreeSIAMGreensFunction<FPType>, 
TList<
      VertexAdd<DDConfSIAM>,
                            TList<VertexRemove<DDConfSIAM>,
                                                            TList<MTMVertexAdd<DDConfSIAM>, NullType>
                                 >
     >
>
{//The init helper for SIAM
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConfSIAM::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConfSIAM>(curparams.U, curparams.beta, 1, alpha);
        moves[1] = new VertexRemove<DDConfSIAM>(curparams.U, curparams.beta, 1, alpha);
	moves[2] = new MTMVertexAdd<DDConfSIAM>(curparams.U, curparams.beta, 1, alpha);
    }
};

template <typename FPType, class DDConfTI>
struct DDQMC_Init_Helper<KondoImpTI<FPType>, TList<VertexAdd<DDConfTI>, TList<VertexRemove<DDConfTI>, NullType> > >
{//The init helper for the Kondo impurity
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConfTI::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConfTI>(curparams.U, curparams.beta, 1, alpha);
        moves[1] = new VertexRemove<DDConfTI>(curparams.U, curparams.beta, 1, alpha);
    }
};

template <typename FPType, class DDConfTI>
struct DDQMC_Init_Helper<KondoImpTI<FPType>, TList<VertexAdd<DDConfTI>, TList<VertexRemove<DDConfTI>, TList<SpinFlip<DDConfTI>, NullType> > > >
{//The init helper for the Kondo impurity with 3 moves...
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConfTI::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConfTI>(curparams.U, curparams.beta, 1, alpha);
        moves[1] = new VertexRemove<DDConfTI>(curparams.U, curparams.beta, 1, alpha);
	moves[2] = new SpinFlip<DDConfTI>(curparams.U, curparams.beta, alpha);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<FreeRealTimeHubbardGreensFunctionLessMem<FPType>,
            TList<ColdAtomVertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> > >
{//the init helper for the ColdAtoms simulations with the Less Mem Greensfunction
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new ColdAtomVertexAdd<DDConf>(curparams.U, 2.0*curparams.t_exp, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.beta, curparams.N, alpha);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<FreeRealTimeHubbardGreensFunction<FPType>,
            TList<ColdAtomVertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> > >
{//The init helper for the ColdAtom simulation with the full tabulating Greensfunction
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        DDQMC_Init_Helper<
        FreeRealTimeHubbardGreensFunctionLessMem<FPType>,
        TList<ColdAtomVertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
        >::init(moves, alpha, curparams);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<FreeRealTimeHubbardGreensFunction<FPType>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> > >
{//The init helper for the Hubbard model, that uses a full tabulating Green's function
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        DDQMC_Init_Helper<
        FreeRealTimeHubbardGreensFunctionLessMem<FPType>,
        TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
        >::init(moves, alpha, curparams);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<
            FreeRealTimeHubbardGreensFunctionLessMem<FPType>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
            >
{//An init helper for the Hubbard model that takes a LessMem Greensfunction
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<
            FreeRealTimeHubbardGreensFunctionLessMem<FPType>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, TList<SpinFlip<DDConf>, NullType> > >
            >
{//An init helper for the Hubbard model that takes a LessMem Greensfunction
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
	moves[2] = new SpinFlip<DDConf>(curparams.U, curparams.contourlen, alpha);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<
            FreeHubbardGreensFunction<FPType>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
            >
{//An init helper for the imaginary time Hubbard model
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<
            RashbaChain<FPType>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
            >
{//An init helper for the imaginary time Hubbard model with Rashba interaction. Might prove useful if one wants to fiddle with other moves.
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<
            ExponentialRashbaChain<FPType>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
            >
{//An init helper for the imaginary time Hubbard model with Rashba interaction and exponentially decaying range. Might prove useful if one wants to fiddle with other moves.
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
    }
};

template <typename FPType, class DDConf>
struct DDQMC_Init_Helper<
            PowerLawRashbaChain<FPType>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
            >
{//An init helper for the imaginary time Hubbard model with Rashba interaction and exponentially decaying range. Might prove useful if one wants to fiddle with other moves.
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
    }
};

template <typename T, bool Has_Spin_Symmetry, class DDConf>
struct DDQMC_Init_Helper<
            ImagModelFromServer<T, Has_Spin_Symmetry>,
            TList<VertexAdd<DDConf>, TList<VertexRemove<DDConf>, NullType> >
            >
{//An init helper for the imaginary time Hubbard models that are supplied by the server
    template <class BasicMove>
    static inline void init(BasicMove *const moves, const Alpha<typename DDConf::FPType>& alpha, const Parameters& curparams)
    {
        moves[0] = new VertexAdd<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
        moves[1] = new VertexRemove<DDConf>(curparams.U, curparams.contourlen, curparams.N, alpha);
    }
};

/**
The Configuration Template for the DDQMC Modell. the Template Parameters are:
_Container the Container in which to store the Vertices(mostly some STL Container like vector or list)
_Vertex The Basic Structure of a Vertex
_GreensFunction The free - particle GreensFunction, that can be evaluated at a vertex
_Observable A Class that can deduce via an evaluate() call the value of the Observable from a given Configuration. It also contains Information about the type of the Observable
_PRNG the container with the Random Number Generator
_MatrixType This specifies which Matrix Library to use. The Interface is expected to match Flos matrix Library
*/
template <template <class tp, class alloc = std::allocator<tp> > class _Container, class _PRNG, typename T, bool SpinSymmetry, class _MatrixType = MTLICCL::Matrix< MTLICCL::/*Dynamic*/MemArray1D< MTLICCL::/*Config_Dynamic*/Config1D<typename ImagModelFromServer<T, SpinSymmetry>::FreeGreensFunctionReturnValueType> > >
>
struct ImagFromServerConfig
{
    typedef ImagModelFromServer<T, SpinSymmetry> GreensFunction;///< the free Greensfunction that we want to use
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType DeterminantType;///< the type of the determinants that occur throughout the program
    typedef typename GreensFunction::Vertex Vertex;///< the Type of the Vertex that we use
    typedef Configuration_Storage<_Container<Vertex>, DeterminantType, double, _MatrixType, GreensFunction::has_Spin_symmetry> Container;///< the container that is used for storing the Vertices
    typedef _PRNG Prng;///< the random number generator structure that we are using now
    typedef _MatrixType Matrix;///< the Type of the Matrix Datatype that is used throughout the program
    typedef double FPType;///<the Floating Point data type that is used throughout the program
};


template <class App>
void Zombie<App>::init()
{
    cout<<"Infecting host..."<<endl;
    typedef double FPType;
    Models model = static_cast<Models>(this->template recvblock<int32_t>());
    bool needsRealtime = (this->template recvblock<char>() != 0);
    //    typedef Prng<MersenneTwister, ToFloat<MersenneTwister, 4294967295U> > prng;
    typedef Prng<Ziff, ToFloat<Ziff, FPType, 2147483648U> > prng;
    //typedef Prng<Coveyou<2147483647>, ToFloat<Marsaglia_Zaman, 4294967295ul> > prng;
    //Configure Modell

    typedef DDQMC_Config<std::vector, prng, FreeRealTimeHubbardGreensFunctionLessMem, FPType> DDConfHubbardChainRealTime;

//    typedef TList<VertexAdd<DDConfHubbardChainRealTime>, TList<VertexRemove<DDConfHubbardChainRealTime>, NullType> > Basic_Moves_Hubbard;
typedef TList<VertexAdd<DDConfHubbardChainRealTime>, TList<VertexRemove<DDConfHubbardChainRealTime>, TList<SpinFlip<DDConfHubbardChainRealTime>, NullType> > > Basic_Moves_Hubbard;
    typedef DDQMC<DDConfHubbardChainRealTime, Basic_Moves_Hubbard> DDQMChubbardchain;

    typedef TList<ColdAtomVertexAdd<DDConfHubbardChainRealTime>, TList<VertexRemove<DDConfHubbardChainRealTime>, NullType> > Basic_Moves_Cold_Atoms;
    typedef DDQMC<DDConfHubbardChainRealTime, Basic_Moves_Cold_Atoms> DDQMC_Cold_Atoms;

    typedef DDQMC_Config<std::vector, prng, FreeSIAMGreensFunction, FPType>  DDConfSIAM;
    typedef TList<VertexAdd<DDConfSIAM>, TList<VertexRemove<DDConfSIAM>, TList<MTMVertexAdd<DDConfSIAM>, NullType> > > Basic_Moves_SIAM;
    typedef DDQMC<DDConfSIAM, Basic_Moves_SIAM> DDQMCSIAM;

    typedef DDQMC_Config<std::vector, prng, FreeHubbardGreensFunction, FPType> DDConfHubbardChainImagTime;

    typedef TList<VertexAdd<DDConfHubbardChainImagTime>, TList<VertexRemove<DDConfHubbardChainImagTime>, NullType> > Basic_Moves_Hubbard_Imag;
    typedef DDQMC<DDConfHubbardChainImagTime, Basic_Moves_Hubbard_Imag> DDQMChubbardchainImag;
    

    typedef ImagFromServerConfig<std::vector, prng, std::complex<FPType>, false> DDConfImagTimeModelFromServerComplexNoSpinSymmetry;
    typedef ImagFromServerConfig<std::vector, prng, std::complex<FPType>, true> DDConfImagTimeModelFromServerComplexSpinSymmetry;
    typedef ImagFromServerConfig<std::vector, prng, FPType, false> DDConfImagTimeModelFromServerRealNoSpinSymmetry;
    typedef ImagFromServerConfig<std::vector, prng, FPType, true> DDConfImagTimeModelFromServerRealSpinSymmetry;
    
    typedef TList<VertexAdd<DDConfImagTimeModelFromServerComplexNoSpinSymmetry>, TList<VertexRemove<DDConfImagTimeModelFromServerComplexNoSpinSymmetry>, NullType> > Moves_ServerModelComplexNoSpinSymmetry;
    typedef TList<VertexAdd<DDConfImagTimeModelFromServerComplexSpinSymmetry>, TList<VertexRemove<DDConfImagTimeModelFromServerComplexSpinSymmetry>, NullType> > Moves_ServerModelComplexSpinSymmetry;
    typedef TList<VertexAdd<DDConfImagTimeModelFromServerRealNoSpinSymmetry>, TList<VertexRemove<DDConfImagTimeModelFromServerRealNoSpinSymmetry>, NullType> > Moves_ServerModelRealNoSpinSymmetry;
    typedef TList<VertexAdd<DDConfImagTimeModelFromServerRealSpinSymmetry>, TList<VertexRemove<DDConfImagTimeModelFromServerRealSpinSymmetry>, NullType> > Moves_ServerModelRealSpinSymmetry;
    
    typedef DDQMC<DDConfImagTimeModelFromServerComplexNoSpinSymmetry, Moves_ServerModelComplexNoSpinSymmetry> DDQMCServerModelComplexNoSpinSymmetry;
    typedef DDQMC<DDConfImagTimeModelFromServerComplexSpinSymmetry, Moves_ServerModelComplexSpinSymmetry> DDQMCServerModelComplexSpinSymmetry;
    typedef DDQMC<DDConfImagTimeModelFromServerRealNoSpinSymmetry, Moves_ServerModelRealNoSpinSymmetry> DDQMCServerModelRealNoSpinSymmetry;
    typedef DDQMC<DDConfImagTimeModelFromServerRealSpinSymmetry, Moves_ServerModelRealSpinSymmetry> DDQMCServerModelRealSpinSymmetry;
    
    typedef DDQMC_Config<std::vector, prng, KondoImpTI, FPType>  DDConfTI;
    typedef TList<VertexAdd<DDConfTI>, TList<VertexRemove<DDConfTI>, NullType> > Basic_Moves_TI;
    typedef DDQMC<DDConfTI, Basic_Moves_TI> DDQMCTI;
    
    typedef DDQMC_Config<std::vector, prng, RashbaChain, FPType>  DDConfRashba;
    typedef TList<VertexAdd<DDConfRashba>, TList<VertexRemove<DDConfRashba>, NullType> > Basic_Moves_Rashba;
    typedef DDQMC<DDConfRashba, Basic_Moves_Rashba> DDQMCRashba;
    
    typedef DDQMC_Config<std::vector, prng, PowerLawRashbaChain, FPType>  DDConfPowerLawRashba;
    typedef TList<VertexAdd<DDConfPowerLawRashba>, TList<VertexRemove<DDConfPowerLawRashba>, NullType> > Basic_Moves_PowerLawRashba;
    typedef DDQMC<DDConfPowerLawRashba, Basic_Moves_PowerLawRashba> DDQMCPowerLawRashba;
    
    typedef DDQMC_Config<std::vector, prng, ExponentialRashbaChain, FPType>  DDConfExponentialRashba;
    typedef TList<VertexAdd<DDConfExponentialRashba>, TList<VertexRemove<DDConfExponentialRashba>, NullType> > Basic_Moves_ExponentialRashba;
    typedef DDQMC<DDConfExponentialRashba, Basic_Moves_ExponentialRashba> DDQMCExponentialRashba;
    
    
    std::cout<<"model: ";
    if (needsRealtime)
    {//here are the models for which we support realtime evolution
        switch (model)
        {
        case HUBBARD_CHAIN:
            std::cout<<"1D Hubbard model in real time"<<std::endl;
//            app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMChubbardchain, Zombie<App> > >(*this, model);
            break;
        case COLD_ATOMS:
            std::cout<<"1D Chain of ultracold atoms"<<std::endl;
//            app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMC_Cold_Atoms, Zombie<App> > >(*this, model);
            break;
        default:
            std::cout<<"Model not recognized! "<<std::endl;
            throw(std::runtime_error("Model not recognized!"));
            break;
        }
    }
    else
    {//here are the models for which we can do imaginary time simulations
        switch (model)
        {
        case HUBBARD_CHAIN:
            std::cout<<"1D Hubbard model in imaginary time"<<std::endl;
            app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMChubbardchainImag, Zombie<App> > >(*this, model);
            break;
        case SIAM:
            std::cout<<"SIAM in imaginary time"<<std::endl;
            app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCSIAM, Zombie<App> > >(*this, model);
            break;
	case IMAG_MODEL_FROM_FILE:
	{
	  std::cout<<"Imaginary time Simulation with free Greensfunction provided by server."<<std::endl;
	  bool has_Spin_symmetry = (this->template recvblock<char>() != 0);
	  bool complexdata = (this->template recvblock<char>() != 0);
	  std::cout<<"Complex Data: "<<(complexdata?"yes":"no")<<std::endl;
	  std::cout<<"Spin Symmetry: "<<(has_Spin_symmetry?"yes":"no")<<std::endl;
	  uint test = has_Spin_symmetry? 1u :0u;
	  test += (complexdata ? 2u: 0u);
	  //first bit is whether the model has spin-symmetry, second bit denotes complex data
	  switch(test)
	  {
 	    case 0:
// 	      app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCServerModelRealNoSpinSymmetry, Zombie<App> > >(*this, model);
 	      break;
	    case 1:
//	      app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCServerModelRealSpinSymmetry, Zombie<App> > >(*this, model);
	      break;
 	    case 2:
// 	      app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCServerModelComplexNoSpinSymmetry, Zombie<App> > >(*this, model);
 	      break;
	    case 3:
//	      app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCServerModelComplexSpinSymmetry, Zombie<App> > >(*this, model);
	      break;
	  }
	  break;
	}
	    case KONDO_IMP_TI:
	      std::cout<<"Anderson impurity in a topological insulator bath in imaginary time"<<std::endl;
//              app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCTI, Zombie<App> > >(*this, model);
            break;
	    case RASHBA_CHAIN:
	      std::cout<<"1D Hubbard chain with Rashba interaction in imaginary time"<<std::endl;
//	      app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCRashba, Zombie<App> > >(*this, model);
	      break;
	   case RASHBA_CHAIN_EXPONENTIAL:
	      std::cout<<"1D Hubbard chain with Rashba interaction in imaginary time and exponential decay"<<std::endl;
//	      app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCExponentialRashba , Zombie<App> > >(*this, model);
	      break;
	   case RASHBA_CHAIN_POWER_LAW:
	      std::cout<<"1D Hubbard chain with Rashba interaction in imaginary time and power law decay"<<std::endl;
//	      app = new QMC_Frame< QMC_Config<prng, Metropolis, DDQMCPowerLawRashba , Zombie<App> > >(*this, model);
	      break;
        default:
            std::cout<<"Model not recognized! "<<std::endl;
            throw(std::runtime_error("Model not recognized!"));
            break;
        }
    }
}

template <class App>
void Zombie<App>::docycle()
{
    app->docycle();
}

template <class App>
void Zombie<App>::run()
{
    unsigned int cnt = 0;
    while (state != CS_TERMINATE)
    {
        if (isdataready())//heck... we received a command!
        {
            char command = recv<char>();
            cout<<"command: "<<(int)command<<endl;
//evaluate command
            switch (command)
            {
            case TERMINATE:
                state = CS_TERMINATE;
                break;
            case INIT:
                delete app;
		sendtoserver<char>(INIT_RECEIVED);
                init();
                break;
            case STOP:
                state = CS_NONE;
		sendtoserver<char>(STOP_RECEIVED);
                break;
            case START:
                state = CS_RUNNING;
                break;
            case RECONNECT:
                cout<<"Reconnect not implemented!"<<endl;
            default:
                cout<<"Command not understood! not changing state"<<endl;
                break;
            };
        }
//act according to state
        switch (state)
        {
        case CS_NONE:
	    std::cout<<"idling!"<<std::endl;
            sleep(1);//sleep a second maybe sth. has happened then.... zzzz
            break;
        case CS_RUNNING:
            docycle();
            cnt++;
            break;
        case CS_TERMINATE:
            //we want to terminate... maybe some cleanup work?
            delete app;
            app = NULL;
            break;
        case CS_RECONNECT:
            cout<<"Reconnect not supported"<<endl;
        default:
            cout<<"I'm in an unknown state!!! :"<<state<<" Aborting!"<<endl;
            exit(-1);
            break;
        };
    }
    return;
}

template <class App>
Zombie<App>::Zombie(int argc, char *argv[]) : ParallelClient(argc, argv), state(CS_NONE), app(NULL)
{
    cout<<"What's that?! Not dead?! My Soul?! My Body?! I'm undead!  AAAAAAAAAAAAAAHHHHHHHHHHH!!!!!   *madness*"<<endl;
//we seem to have established a connection
    return;
}

int main(int argc, char *argv[])
{
    cout<<"Compiled on "<<__DATE__<<" at "<<__TIME__<<'\n';
    cout<<"Finally ...... eternal peace has arrived......."<<endl;
    Zombie<QMC_Frame_Base> zombie(argc, argv);
    zombie.run();
    return 0;
}
