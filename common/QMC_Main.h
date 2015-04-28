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
#include <stdint.h>
#include <queue>
#include "Zombie.h"
#include "ddqmc.h"
#include "ObservableContainer.h"
#include "AverageSign.h"
#include "observables.h"
#include "Parameters.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**
The Monte Carlo class that encapsulates all Monte Carlo stuff. It can generate a new configuration
and you can retrieve a configuration
*/
template <class Sim, class Configuration, class Sampling_Algo, class PRNG>
class MonteCarlo
{
public:
    /**
    This function tells the Monte-Carlo class to perform one update move
    */
    void generateConfiguration();
    /**
    This class retrieves a configuration
    @return the new configuration
    */
    const Configuration& getConfiguration();
    /**
    This constructs the Monte Carlo class. The parameters are read from the Parameters object
    @param params the object with the parameters
    */
    MonteCarlo ( const Parameters& params );
    ~MonteCarlo() {std::cout<<"Move 0 accepted: "<<move_accepted[0]/static_cast<double>(nr_of_generated_confs)<<"   Move 1 accepted: "<<move_accepted[1]/static_cast<double>(nr_of_generated_confs)<<std::endl;}
private:
    PRNG prng;///< our local RNG
public:
    const Alpha<typename Sim::FPType> alpha;///< the Alpha Object that is used throughout the project
private:
    Sim sim;///< the concrete Monte-Carlo implementation
    Configuration config;///< the local storage for a configuration
    MonteCarlo(const MonteCarlo&);
    MonteCarlo& operator=(const MonteCarlo&);
    long long int nr_of_generated_confs;
    long long int move_accepted[2];
    int wait;
    int current_wait;
};

template <class Sim, class Configuration, class Sampling_Algo, class PRNG>
MonteCarlo<Sim, Configuration, Sampling_Algo, PRNG>::MonteCarlo ( const Parameters& params ) : prng ( rand(), rand() ), //seed our PRNGs with the system PRNG
        alpha ( params.delta ), sim ( alpha, params ), config ( 1.0 ), wait(1), current_wait(0)
{
    /*    Configuration temp;
        typedef typename Configuration::value_type Vertex;
        temp.push_back(Vertex(0.95131, 0, 1));
        temp.push_back(Vertex(3.42178, 0, 1));
        temp.push_back(Vertex(0.139423, 1, -1));
        temp.push_back(Vertex(0.809612, 1, -1));
        temp.push_back(Vertex(2.13979, 0, 1));
        temp.push_back(Vertex(3.39672, 0, 1));
        temp.push_back(Vertex(4.6518, 1, 1));
        temp.push_back(Vertex(1.02577, 1, -1));
        temp.template reevaluate<typename Sim::Config::GreensFunction> ( alpha );
        cout<<temp.phase<<endl;
        exit(-1);*/
    nr_of_generated_confs = 0;
    move_accepted[0] = 0;
    move_accepted[1] = 0;
    //more debugging cruft.. the phase should be (0.521583,-0.853201) within 3 decimals using time-dependent Hubbard model
}

template <class Sim, class Configuration, class Sampling_Algo, class PRNG>
void MonteCarlo<Sim, Configuration, Sampling_Algo, PRNG>::generateConfiguration()
{
    typedef typename Sim::FPType FPType;
    FPType test = prng.rndfloat();
    unsigned int k = 0;
    //determine which move will be executed
    while ( test > sim.moves[k]->CumulativeRatio ) ++k;
    //Let the Move generate a proposal Configuration
    sim.moves[k]->generateProposalConfiguration ( config );
    //Determine if the Proposal Configurations gets accepted
    FPType t2 = Sampling_Algo::F ( sim.moves[k]->generateProposalPropability() );//get the propability of a move
    FPType t1 = prng.rndfloat();//generate the random number
//    std::cout<<"k = 0 <=> addition k = 1 <=> removal "<<std::endl;
//    std::cout<<"rnd:"<<t1<<" weight: "<<t2<<"  => Move ["<<k<<"] "<<(t1<t2?"accepted":"rejected")<<std::endl;
    if ( t1 < t2 )
    {
      move_accepted[k]++;
      ++current_wait;
        config = sim.moves[k]->returnProposalConfiguration( config );//If yes, copy the Configuration
    }
//   std::cout<<current_wait<<" "<<wait<<std::endl;
//     if((current_wait > wait) && config.size() > 0)
//     {
//       //do one MTM step
//       sim.moves[2]->generateProposalConfiguration(config);
//       FPType t2 = Sampling_Algo::F ( sim.moves[2]->generateProposalPropability() );//get the propability of the MTM Addition
//       bool move_accepted = false;
//       if (t2 >= 1)
// 	move_accepted = true;
//       else
//       {
// 	if (prng.rndfloat() < t2 )
// 	  move_accepted = true;
//       }
//       
//       if (move_accepted)
// 	config = sim.moves[2]->returnProposalConfiguration( config );//If yes, copy the Configuration
//       wait = config.size();
//       current_wait = 0;
// //     std::cout<<wait<<std::endl;
//      //exit(-1);
//     }
/*    std::cout<<config<<std::endl;
    std::cout<<"Inverse: "<<std::endl;
    std::cout<<config.matcont.mat<<std::endl;
    std::cout<<"matrix."<<std::endl;
    std::cout<<MTLICCL::inverse(config.matcont.mat)<<std::endl;
    std::cout<<"Now reevaluating"<<std::endl;
    config.template reevaluate<typename Sim::Config::GreensFunction>(alpha);
    std::cout<<"Inverse: "<<std::endl;
    std::cout<<config.matcont.mat<<std::endl;
    std::cout<<"matrix."<<std::endl;
    std::cout<<MTLICCL::inverse(config.matcont.mat)<<std::endl;
    std::cout<<"Now reevaluating"<<std::endl;*/
++nr_of_generated_confs;
    return;
}

template <class Sim, class Configuration, class Sampling_Algo, class PRNG>
const Configuration& MonteCarlo<Sim, Configuration, Sampling_Algo, PRNG>::getConfiguration()
{
  wait = config.size();
    if (config.template needsreinit<typename Sim::Config::GreensFunction> ( alpha ))
    {
        std::cout<<"rebuilding Matrices!!"<<std::endl;
        config.template reevaluate<typename Sim::Config::GreensFunction> ( alpha );
    }
    return config;
}

class QMC_Frame_Base
{
public:
    virtual void docycle() = 0;
    virtual ~QMC_Frame_Base() {}
};

/**
This Class represents the Basic QMC Framework. It takes a config template that configures it.
It's main purpose is to build the qmc_main function.
*/
template <class Config>
class QMC_Frame : public QMC_Frame_Base
{
public:
    /**
    Basic Monte Carlo loop. Evaluates the model given by Config::App. One warm-up-cycle with N_wait steps is done at the beginning.
    The Waiting Time(how many configurations) between consecutive measurements cycles is retrieved from the server
    */
    typedef typename Config::App::Config Sim_Config;//For easy Access to the Config Repository of the modell
    typedef typename Sim_Config::Container Configuration;//a typedef for the data type of the Configuration
    typedef typename Sim_Config::GreensFunction GreensFunction;// a typedef for the Greensfunction
    typedef typename Sim_Config::GreensFunction::FreeGreensFunctionReturnValueType GFRetVal;//the type of the sign
    typedef QMC_Frame<Config> MyType;
    typedef AverageSign<Configuration, GFRetVal> AverageSignType;// a typedef for the DataType of the AverageSign
    typedef typename Config::Comm Comm;// a typedef for the Communication object
    QMC_Frame (Comm&, Models);
    inline void docycle ();
    inline ~QMC_Frame();
private:
    const Parameters curparams;
    struct SysPrngIniter
    {
        SysPrngIniter ( unsigned int x )
        {
            srand ( x );
        }
    } sysprnginiter;
    ObservableContainer<Configuration, GreensFunction, typename Sim_Config::FPType> observables;///< A cache where all observables are stored
    const unsigned int nr_of_walkers;
    MonteCarlo<typename Config::App, Configuration, typename Config::Sampling_Algo, typename Config::PRNG>** montecarlo;///< A variable that holds the state of a Monte-Carlo walker
    bool* walkerused;
    AverageSignType averagephase;///< the averagesign is something fundamental... thus we store it extra
    std::queue<Configuration> fifo;
    const unsigned int maxlevel;
    const unsigned int minlevel;
    int requestnewWalkerID();
    void freeWalker(int);
#ifdef _OPENMP
    bool generateconfig;
#endif
    Comm& net;
};

template <class Config>
QMC_Frame<Config>::~QMC_Frame()
{
    delete [] walkerused;
    for(unsigned int k = 0; k < nr_of_walkers; ++k )
      delete montecarlo[k];
    delete [] montecarlo;
    Sim_Config::GreensFunction::tidyup();
}

static int omphelper()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

template <class Config>
void QMC_Frame<Config>::freeWalker(int id)
{
    walkerused[id] = false;
}

template <class Config>
int QMC_Frame<Config>::requestnewWalkerID()
{
    int tries = 1000;
    int retval = -1;
    do
    {
        int test = rand() % nr_of_walkers;
        if (walkerused[test] == false)
        {
            retval = test;
            walkerused[test] = true;
        }
        tries--;
    } while (( retval == -1) && (tries > 0));
    return retval;
}

template <class Obs_Config, class GreensFunction, bool has_GIomegan = false>
struct HasGiomegan
{
    template<class Comm, class CFG, class Container>
    static inline bool add_observables(const std::string&, Comm& z, CFG& curparams, Container&) {}
};

template <class Obs_Config, class GreensFunction, bool has_GIomegan = false>
struct HasTRS
{
    template<class Comm, class CFG, class Container>
    static inline bool add_observables(const std::string&, Comm& z, CFG& curparams, Container&) {}
};

template <class Obs_Config, class GreensFunction, bool has_GIomegan = false, bool Is_Impurity = false, bool Spin_Symmetry = false>
struct UseLocalBathGreensfunctions
{
    template<class Comm, class CFG, class Container>
    static inline bool add_observables(const std::string&, Comm& z, CFG& curparams, Container&) {}
};

template <class Obs_Config, class GreensFunction>
struct HasGiomegan<Obs_Config, GreensFunction, true>
{
    template<class Comm, class CFG, class Container>
    static inline bool add_observables(const std::string& obs, Comm& z, CFG& curparams, Container& cont)
    {
      bool retval = true;
        if (obs == "MatsubaraFrequencyGreensfunction")
            cont.add_Observable(new MatsubaraFrequencyGreensfunction<Obs_Config, GreensFunction>(z, curparams));
	else
	if (obs == "MatsubaraFrequencyGreensfunction_averaged")
            cont.add_Observable(new MatsubaraFrequencyGreensfunction<Obs_Config, GreensFunction>(z, curparams));
	else
//        if (obs == "MatsubaraFrequencyGreensfunctionDerivative")
//            cont.add_Observable(new MatsubaraFrequencyGreensfunctionDerivative<Obs_Config, GreensFunction>(z, curparams));
//	  else
//        if (obs == "Conductance")
//            cont.add_Observable(new Conductance<Obs_Config, GreensFunction>(z, curparams));
//	    else
if(obs == "SpinSusceptibility_X")
  cont.add_Observable(new SpinSusceptibility_X<Obs_Config, GreensFunction>(z, curparams));
else if(obs == "SpinSusceptibility_Z")
  cont.add_Observable(new SpinSusceptibility_Z<Obs_Config, GreensFunction>(z, curparams));
else retval = false;
return retval;
    }
};

template <class Obs_Config, class GreensFunction>
struct HasTRS<Obs_Config, GreensFunction, true>
{
    template<class Comm, class CFG, class Container>
    static inline bool add_observables(const std::string& obs, Comm& z, CFG& curparams, Container& cont)
    {
      bool retval = true;
              if (obs == "DiagonalGreensfunction_kspace")
            cont.add_Observable(new DiagonalGreensfunction_kspace<Obs_Config>(z, curparams));
	    else
              if (obs == "OffDiagonalGreensfunction_kspace")
            cont.add_Observable(new OffDiagonalGreensfunction_kspace<Obs_Config>(z, curparams));
	    else
	      retval = false;
	    return retval;
    }
};

template <class Obs_Config, class GreensFunction>
struct UseLocalBathGreensfunctions<Obs_Config, GreensFunction, true, true, false>
{
    template<class Comm, class CFG, class Container>
    static inline bool add_observables(const std::string& obs, Comm& z, CFG& curparams, Container& cont)
    {
      bool retval = true;
      if(obs == "LocalBathGreensfunctions_up")
        cont.add_Observable(new LocalBathGreensfunctions<Obs_Config, GreensFunction, UP>(z, curparams));
      else
      if(obs == "LocalBathGreensfunctions_down")
        cont.add_Observable(new LocalBathGreensfunctions<Obs_Config, GreensFunction, DOWN>(z, curparams));
      else
      if(obs == "LocalBathGreensfunctions_averaged")
        cont.add_Observable(new LocalBathGreensfunctions_averaged<Obs_Config, GreensFunction>(z, curparams));
      else
      if(obs == "KondoCloud_Z")
        cont.add_Observable(new KondoCloud_Z<Obs_Config, GreensFunction, DOWN>(z, curparams));
      else
      if(obs == "KondoCloud_X")
        cont.add_Observable(new KondoCloud_X<Obs_Config, GreensFunction, DOWN>(z, curparams));
      else
      if(obs == "Hybridization_up")
        cont.add_Observable(new Hybridization<Obs_Config, GreensFunction, UP>(z, curparams));
      else
      if(obs == "Hybridization_down")
        cont.add_Observable(new Hybridization<Obs_Config, GreensFunction, DOWN>(z, curparams));
      else
	retval = false;
      return retval;
    }
};

template <class Config>
QMC_Frame<Config>::QMC_Frame (Comm& z, Models model)
        :
        curparams(z, model), //retrieve all parameters from the master
        sysprnginiter ( z.template recvblock<int32_t>() ),//seed the system PRNG
        observables(curparams),
        nr_of_walkers(omphelper()),
        montecarlo (new MonteCarlo<typename Config::App, Configuration, typename Config::Sampling_Algo, typename Config::PRNG>*[nr_of_walkers] ),
        walkerused(new bool[nr_of_walkers]),
        maxlevel(4*curparams.prebins), minlevel(curparams.prebins + 1),
        net(z)
{
    //initialize the Greensfunction tables
    std::cout<<"creating Greensfunction tables"<<std::endl;
    GreensFunction::init(curparams);
    std::cout<<"Generating space for "<< nr_of_walkers << " walkers"<<std::endl;
    for (unsigned int k = 0; k < nr_of_walkers; ++k)
        walkerused[k] = false;
    uint32_t nrof_observables = z.template recvblock<uint32_t>();
    typedef Observable_Parameters<Configuration, GFRetVal, Comm> Obs_Config;
    for ( unsigned int k = 0; k < nrof_observables; ++k )
    {
        const std::string obs ( z.template recvblock<std::string>() );
        if ( obs == "AverageOrder" )
            this->observables.add_Observable ( new AverageOrder<Obs_Config>( z ) );
        else
            if ( obs == "ParticleNumber" )
                this->observables.add_Observable ( new ParticleNumber<Obs_Config>( z, curparams ) );
            else
                if ( obs == "TotalDoubleOccupancy" )
                    this->observables.add_Observable ( new TotalDoubleOccupancy<Obs_Config>( z, curparams ) );
                else
                    if ( obs == "KineticEnergy" )
                        this->observables.add_Observable ( new KineticEnergy<Obs_Config>( z, curparams ) );
                    else
                        if ( obs == "SpinSpinCorrelation" )
                            this->observables.add_Observable ( new SpinSpinCorrelation<Obs_Config>( z, curparams ) );
                        else
                            if ( obs == "ChargeChargeCorrelation" )
                                this->observables.add_Observable ( new ChargeChargeCorrelation<Obs_Config> ( z, curparams ) );
                            else
                                if ( obs == "kSpaceDensity" )
                                    this->observables.add_Observable ( new kSpaceDensity<Obs_Config> ( z, curparams ) );
                                else
                                    if ( obs == "Greensfunction" )
                                        this->observables.add_Observable ( new Greensfunction<Obs_Config>( z, curparams ) );
                                    else
                                        if ( obs == "EtaPairing" )
                                            this->observables.add_Observable ( new EtaPairing<Obs_Config> ( z, curparams ) );
                                        else
                                            if ( obs == "LocalDensityVariance" )
                                                this->observables.add_Observable ( new LocalDensityVariance<Obs_Config>( z, curparams ) );
                                            else
                                                if ( obs == "DensityDensityStructureFactor" )
                                                    this->observables.add_Observable ( new DensityDensityStructureFactor<Obs_Config> ( z, curparams ) );
                                                else
                                                    if ( obs == "DoubleOccupancy" )
                                                        this->observables.add_Observable ( new DoubleOccupancy<Obs_Config>( z, curparams ) );
                                                    else
                                                        if ( obs == "ChargeChargeCorrelatedPart" )
                                                            this->observables.add_Observable ( new ChargeChargeCorrelatedPart<Obs_Config> ( z, curparams ) );
                                                        else
                                                            if ( obs == "SpinSpinCorrelatedPart" )
                                                                this->observables.add_Observable ( new SpinSpinCorrelatedPart<Obs_Config> ( z, curparams ) );
                                                            else
                                                                if ( obs == "EtaPairing_Real" )
                                                                    this->observables.add_Observable ( new EtaPairing_Real<Obs_Config> ( z, curparams ) );
                                                                else
                                                                    if ( obs == "Magnetization" )
                                                                        this->observables.add_Observable ( new Magnetization<Obs_Config> ( z, curparams ) );
                                                                    else
								      if ( obs == "SpinSpinCorrelatedPart_Y" )
                                                                        this->observables.add_Observable ( new SpinSpinCorrelatedPart_Y<Obs_Config> ( z, curparams ) );
								      else
								        if ( obs == "SpinSpinCorrelatedPart_X" )
                                                                          this->observables.add_Observable ( new SpinSpinCorrelatedPart_X<Obs_Config> ( z, curparams ) );
								        else
                                                                    if(curparams.t_exp == 0.0)
                                                                    {
								      if ( obs == "ImaginarySpinSpinCorrelation_Z" )
                                                                          this->observables.add_Observable ( new ImaginarySpinSpinCorrelation_Z<Obs_Config> ( z, curparams ) );
								        else
									  if (obs == "SpinSpinCorrelation_kspace")
									    this->observables.add_Observable( new SpinChargeParent<Obs_Config, -1> ( z, curparams ) );
									  else 
									    if (obs == "ChargeChargeCorrelation_kspace")
									      this->observables.add_Observable( new SpinChargeParent<Obs_Config, 1> ( z, curparams ) );
									    else
									    	if (obs == "SplusSminus_kspace")
									      this->observables.add_Observable( new SplusSminus<Obs_Config> ( z, curparams ) );
									    else
									      if (obs == "Gplusplus_kspace")
									      this->observables.add_Observable( new TRI_Greensfunction_kspace<Obs_Config, 1> ( z, curparams ) );
									    else
									      if (obs == "Gminusminus_kspace")
									      this->observables.add_Observable( new TRI_Greensfunction_kspace<Obs_Config, -1> ( z, curparams ) );
									    else
                                                                             if (obs == "SmoothImaginaryGreensfunction")
                                                                                this->observables.add_Observable ( new SmoothImaginaryGreensfunction<Obs_Config>( z, curparams ) );
                                                                             else
									     if (obs == "SmoothImaginaryGreensfunction_averaged")
                                                                                this->observables.add_Observable ( new SmoothImaginaryGreensfunction_averaged<Obs_Config>( z, curparams ) );
									     else
									     {
									     bool found = HasGiomegan< Obs_Config, GreensFunction, GreensFunction::has_Giomega>::template add_observables(obs, z, curparams, this->observables) ||
									     HasTRS< Obs_Config, GreensFunction, GreensFunction::has_TRS>::template add_observables(obs, z, curparams, this->observables) ||
									     UseLocalBathGreensfunctions<Obs_Config, GreensFunction, GreensFunction::has_Giomega,
									     GreensFunction::Is_Impurity_model, GreensFunction::has_Spin_symmetry
									     >::template add_observables(obs, z, curparams, this->observables);
									     if(!found)
                                                                              {
                                                                            std::cout<<"Observable "<<obs<<" not known! Aborting"<<std::endl;
                                                                            exit ( -1 );
                                                                        }
									     }
                                                                    }
                                                                        else
                                                                        {
                                                                            std::cout<<"Observable "<<obs<<" not known! Aborting"<<std::endl;
                                                                            exit ( -1 );
                                                                        }
    }
    observables.dryrun();
    std::cout<<"Nr of Entries that will be cached: "<<observables.num_Cache_Entries() <<std::endl;
    //observables.dump_cache();
    //let's create the walkers
    for (unsigned int k = 0; k < nr_of_walkers; ++k)
    {
        std::cout<<"Generating Walker "<<k<<std::endl;
        montecarlo[k] = new MonteCarlo<typename Config::App, Configuration, typename Config::Sampling_Algo, typename Config::PRNG>(curparams);
    }
    //the warm up Loop
    std::cout<<"Warming Up"<<std::endl;
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < static_cast<int>(nr_of_walkers); ++k)
    {
        for ( unsigned int cnt = 0; cnt < /*20 * curparams.N_wait*/2000; ++cnt )
	{
            montecarlo[k]->generateConfiguration();
	}
    }
//make sure that the output is not messed up
    for (unsigned int k = 0; k < nr_of_walkers; ++k)
    {
        std::cout<<"Warmup - Phase finished... Configuration of walker "<<k<<" is now:"<<'\n';
        auto res = montecarlo[k]->getConfiguration();
        std::cout<<res<<'\n';
// //        std::cout<<res.matcont.up<<std::endl;
// //	std::cout<<(montecarlo[k]->getConfiguration().matcont.needsreinit< Alpha>()?"yes":"no")<<std::endl;
// 	std::cout<<det(res.matcont.up)<<std::endl;
// 	std::cout<<det(res.matcont.down)<<std::endl;
// 	std::cout<<res.matcont.upratio<<std::endl;
// 	std::cout<<res.matcont.downratio<<std::endl;
// 	std::cout<<res.matcont.getRatio()<<std::endl;
// 	res.template reevaluate<GreensFunction> (montecarlo[k]->alpha);
// //	        std::cout<<res.matcont.up<<std::endl;
// //	std::cout<<(montecarlo[k]->getConfiguration().matcont.needsreinit< Alpha>()?"yes":"no")<<std::endl;
// 	std::cout<<det(res.matcont.up)<<std::endl;
// 	std::cout<<det(res.matcont.down)<<std::endl;
// 	std::cout<<res.matcont.upratio<<std::endl;
// 	std::cout<<res.matcont.downratio<<std::endl;
// 	std::cout<<res.matcont.getRatio()<<std::endl;
    }
#ifdef _OPENMP
    omp_set_nested ( true );
    if ( !omp_get_nested() ) std::cout<<"It seems OpenMP nesting is not supported!"<<std::endl;
    generateconfig = true;
#endif
    return;
}

template <class Config>
void QMC_Frame<Config>::docycle ()
{
//   memctr = 0;
//   memsize = 0;
//   start = clock();
   const uint32_t& prebins = curparams.prebins;
#ifdef _OPENMP
// #pragma omp parallel
//     {
// #pragma omp single
//         {//specify that only the master thread generates configs
//             while (generateconfig == true)
//             {
//                 int w = requestnewWalkerID();//master tries to get a new ID for a walker
//                 if (w != -1)//a walker is available, so let's generate a task that uses it
//                 {
//                     std::cout<<"generating thread for walker: "<<w<<std::endl;
// #pragma omp task firstprivate(w)
//                     {//begin of code of task
//                         std::cout<<"running  walker "<<w <<"!"<<std::endl;
// 			MonteCarlo<typename Config::App, Configuration, typename Config::Sampling_Algo, typename Config::PRNG>& mymc(*(montecarlo[w]));
//                         for (unsigned int k = 0; k < curparams.N_wait/prebins; ++k)
//                             mymc.generateConfiguration();
// #pragma omp critical
//                         {
// //			    std::cout<<"here walker "<<w<<std::endl;
//                             fifo.push(mymc.getConfiguration());
// //			    std::cout<<"stored config:"<<std::endl;
// //                            std::cout<<fifo.back()<<std::endl;
//                         }
//                         freeWalker(w);
//                         std::cout<<"walker "<<w<<" freed"<<std::endl;
//                     }//end of code of task
//                 }
//                 else
//                     std::cout<<"No walker available!"<<std::endl;//the master thread wastes time in here....
#pragma omp parallel for
for(int w = 0; w < nr_of_walkers; ++w)
{
        for (unsigned int k = 0; k < curparams.N_wait/prebins; ++k)
             montecarlo[w]->generateConfiguration();
}
for(int w = 0; w < nr_of_walkers; ++w)
{
  fifo.push(montecarlo[w]->getConfiguration());
}

                if (fifo.size() > maxlevel)
                    generateconfig = false;
//                if (w == -1)
//                    generateconfig = false;
//            }

            //master checks whether the queue is full
            if (fifo.size() >= maxlevel)
            {
                std::cout<<"maximum number of configs reached!"<<std::endl;
                generateconfig = false;//we reached the maximum number of configurations that we want to cache, so let's stop further creation of tasks
            }
 //           }//end of master region
//    }//end of the parallel region
//now enter measurement part
            if (fifo.size() > minlevel)
            {
	      while(fifo.size() > minlevel)
	      {
                for (unsigned int i = 0; i < prebins; ++i)
                {
                    Configuration config;
#pragma omp critical
                    {//protect access to the queue
                        std::cout<<"configurations left: "<<fifo.size()<<std::endl;
                        config = fifo.front();
                        fifo.pop();
                    }
                    //measure the average sign
                    averagephase.evaluate(config);
                    //Evaluate the Value of the observables for the Configuration
                    observables.measure(config);
                }
                averagephase.template sendData<Comm>(net);
                observables.executefunction(&ObservableBase<Configuration, GFRetVal>::senddata);
                generateconfig = true;
	      }
            }
            else
            {//we used most of the available configurations
                std::cout<<"minimum level reached!"<<std::endl;
                generateconfig = true;
            }
//        }//end of master region
//    }//end of the parallel region
#else
    for (unsigned int j = 0; j < prebins; ++j)
    {
        //now follows the real Loop that generates configurations. we propose N_wait/prebins configurations between every measurement
        for (unsigned int cnt = 0; cnt < curparams.N_wait/prebins; ++cnt)
            montecarlo[0]->generateConfiguration();
        //reevaluate the matrices of the configuration right before doing measurements. This is done by the MonteCarlo Class
        Configuration configuration(montecarlo[0]->getConfiguration());
        std::cout<<"[QMC_Main] Configuration_Length: "<<configuration.size()<<std::endl;
        std::cout<<"[QMC_Main] Configuration phase: "<<configuration.phase<<std::endl;
        //measure the average sign
        averagephase.evaluate(configuration);
        //Evaluate the Value of the observables for the Configuration
        observables.measure(configuration);
    }
//    std::cout<<"sending data to server"<<std::endl;
    averagephase.template sendData<Comm>(net);
    observables.executefunction(&ObservableBase<Configuration, GFRetVal>::senddata);
#endif
//     endtime = clock();
//     std::cout<<"Allocations: "<<memctr<<std::endl;
//     std::cout<<"Average mem: "<<static_cast<float>(memsize)/static_cast<float>(memctr)/1024.0<<" kB"<<std::endl;
//     std::cout<<"Runtime: "<<((endtime - start) * 1000)/CLOCKS_PER_SEC<<" ms"<<std::endl;
//     exit(-1);
}

/**
This is the Configuration Structure for the QMC_Main Object. The template parameters are as follows:
     _PRNG The PseudoRandomNumberGenerator to use. It must have the methods rndfloat() and rndint()
     _Sampling_Algo The Sampling Algorithm to use for Importance Sampling. Must have a method T F(T)
     _App the actual Model to simulate
*/
template <class T_PRNG, class T_Sampling_Algo, class T_App, class Network>
struct QMC_Config
{
    typedef T_PRNG PRNG;///< A typedef for the RNGs that is used throughout the Program
    typedef T_Sampling_Algo  Sampling_Algo;///< a typedef for the sampling Algorithm that we use
    typedef T_App App;///< A typedef for the actual model that we evaluate
    typedef Network Comm;///< this is the type of the object where we retrieve data from
};
