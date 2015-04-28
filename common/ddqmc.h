#ifndef DDQMC_H
#define DDQMC_H
#include <iostream>
#include <complex>
#include "generalPhysics.h"
#include "Vertex.h"
#include "Parameters.h"
#include "TList.h"
#include "Configuration_Storage.h"

#include <cassert>

/**
A helper function to update the members of a partially updated configuration
This is the spin symmetric version.
@param configuration the configuration to update
*/
template<typename DetType, typename FPType>
static inline void updateWeights(DetType& ratioA, FPType& interpretableratio, DetType& incphase, std::pair<DetType, DetType> ratio)
{
    ratioA = ratio.first * ratio.second;
    interpretableratio = std::abs(ratioA);//generate the weight that we can interpret as propability. it's independent of all the phase-factors
    incphase *= ratioA/interpretableratio;
    return;
}

/**
A helper function to update the members of a partially updated configuration.
This Version is called in the general case
@param configuration the configuration to update
*/
template<typename DetType, typename FPType, typename T>
static inline void updateWeights(DetType& ratioA, FPType& interpretableratio, DetType& incphase, const T& ratio)
{
    ratioA = MTLICCL::det(ratio);
//    std::cout<<"[updateWeights] ratioA = "<<ratioA<<std::endl;
    interpretableratio = std::abs(ratioA);//generate the weight that we can interpret as propability. it's independent of all the phase-factors
    incphase *= ratioA/interpretableratio;
    return;
}

/**
The Configuration Template for the DDQMC Modell. the Template Parameters are:
_Container the Container in which to store the Vertices(mostly some STL Container like vector or list)
_Vertex The Basic Structure of a Vertex
_GreensFunction The free - particle GreensFunction, that can be evaluated at a vertex
_Observable A Class that can deduce via an evaluate() call the value of the Observable from a given Configuration. It also contains Information about the type of the Observable
_PRNG the container with the Random Number Generator
_MatrixType This specifies which Matrix Library to use. The Interface is expected to match Flos matrix Library
*/
template <template <class tp, class alloc = std::allocator<tp> > class _Container, class _PRNG, template<class> class _GreensFunction, typename _FPType = float, class _MatrixType = MTLICCL::Matrix< MTLICCL::/*Dynamic*/MemArray1D< MTLICCL::/*Config_Dynamic*/Config1D<typename _GreensFunction<_FPType>::FreeGreensFunctionReturnValueType> > >
>
struct DDQMC_Config
{
    typedef typename _GreensFunction< _FPType>::FreeGreensFunctionReturnValueType DeterminantType;///< the type of the determinants that occur throughout the program
    typedef typename _GreensFunction< _FPType>::Vertex Vertex;///< the Type of the Vertex that we use
    typedef _GreensFunction<_FPType> GreensFunction;///< the free Greensfunction that we want to use
    typedef Configuration_Storage<_Container<Vertex>, DeterminantType, _FPType, _MatrixType, GreensFunction::has_Spin_symmetry> Container;///< the container that is used for storing the Vertices
    typedef _PRNG Prng;///< the random number generator structure that we are using now
    typedef _MatrixType Matrix;///< the Type of the Matrix Datatype that is used throughout the program
    typedef _FPType FPType;///<the Floating Point data type that is used throughout the program
};

/**
Base Class for a move. generateProposalConfiguration(), generateProposalPropability(), returnProposalConfiguration() are the basic interface to the moves. As every move has a propability the Cumulative Ratio of the Moves gets also stored here. Via _Config the Configuration of the model parameters is passed to the moves.
*/
template <class _Config>
class GeneralMove
{
public:
    typedef _Config Config;///< a typedef to access the configuration
    typedef typename Config::Container Configuration;///<a typedef for the used container
    typedef typename _Config::Prng Prng;///< a typedef for the used RNG
    /**
    This function is supposed to generate the changed configuration
    @param c the old configuration
    */
    virtual void generateProposalConfiguration(const Configuration& c) = 0;
    /**
    This function generates a propability that is plugged into the sampling algorithm
    @return the propability of the proposed move
    */
    virtual typename Config::FPType generateProposalPropability() = 0;
    /**
    If the new configuration is accepted, this function returns it.
    @param c the old configuration
    @return the configuration that is proposed
    */
    virtual Configuration& returnProposalConfiguration(const Configuration& c) MTL_CONST_FUNCTION;
    /**
    the destructor of the moves
    */
    virtual ~GeneralMove() {}
    /**
    the constructor that any derived constructor should call
    @param cr the cumulative Ratio for the move
    @param a the alpha object that the move uses
    */
    inline GeneralMove(float cr, const Alpha<typename Config::FPType>& a);
    const float CumulativeRatio;///< The CumulativeRatio of every Move
protected:
    typename _Config::DeterminantType ratio;///< here we store the cumulative ratio of the moves
    typename _Config::FPType interpretableratio;
    Configuration proposalconfiguration;///< a place to store the configuration that will be proposed
    const Alpha<typename Config::FPType>& alpha;///< the alpha Object that the move uses
private:
};

template<class _Config>
GeneralMove<_Config>::GeneralMove(float cr, const Alpha<typename Config::FPType>& a) : CumulativeRatio(cr), alpha(a)
{
    return;
}

template<class _Config>
typename _Config::Container& GeneralMove<_Config>::returnProposalConfiguration(const Configuration& c)
{
    return this->proposalconfiguration;
}

/**
The sole purpose of this class is to determine the sign of the configuration.
The integer template parameter is set to zero. This signifies that we request an imaginary time QMC.
*/
template<typename GreensFunction, class Configuration>
struct PhaseEvaluator<GreensFunction, Configuration, 0>
{
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType RetType;
    static inline void determinePhase(Configuration& configuration) throw()
    {
        configuration.phase = configuration.incphase;
        if ((configuration.size() & 1)) //add the sign of the weight for the Hubbard-Modell. This is determined by a test for odd or even
            configuration.phase *= static_cast<RetType>(-1);
    }
};

/**
The sole purpose of this class is to determine the sign of the configuration.
The integer template parameter is set to one. This signifies that we request a realtime handling of the sign.
*/
template<typename GreensFunction, class Configuration>
struct PhaseEvaluator<GreensFunction, Configuration, 1>
{
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType RetType;
    static inline void determinePhase(Configuration& c) throw()
    {
        RetType retval;
        switch (c.size() % 4)//determine the phases which are (-i)^n
        {
        case 0:
            retval = 1.0;
            break;
        case 1:
            retval = RetType(0.0,-1.0);
            break;
        case 2:
            retval = -1.0;
            break;
        case 3:
            retval = RetType(0.0,1.0);
            break;
        };
        for (typename Configuration::const_iterator it = c.begin(); it != c.end(); ++it)
        {
            if (it->tau >= GreensFunction::t_exp)
            {
                if (it->tau >= 2.0 * GreensFunction::t_exp)
                    retval *= RetType(0.0,-1.0);//s is on imaginary axis
                else
                    retval *= -1.0;//s is on backward branch
            }
        }
        c.phase = retval * c.incphase;//the total phase are the total incremental phases times the various factors
    }
};

/**
A class that encapsulates the steps for adding a Vertex in the ColdAtomscase
*/
template<class Config>
class ColdAtomVertexAdd : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< a typedef to access the Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< a typedef for the used PRNG
    typedef typename Config::Vertex VA_Vertex;///< the used Vertex Type
    typedef typename Config::FPType FPType;///< the used floating point type
    typedef typename Config::Matrix MatType;
    typedef typename Config::GreensFunction::FreeGreensFunctionReturnValueType GFType;
public:
    /**
    This is the Constructor for the move to add Vertices.
    @param _u the Hubbard U
    @param _N the Number of sites
    @param a the used alpha object
    */
    ColdAtomVertexAdd(FPType _u, FPType low, FPType max, unsigned int _N, const Alpha<FPType>& a) : GeneralMove<Config>(0.5/*the cumulative ratio of this move*/, a), u(_u),lowerpos(low), upperpos(max), contourpartlen(max-low), N(_N), prng(rand(), rand()), prefactor(u * contourpartlen * static_cast<FPType>(N))
    {}
    ~ColdAtomVertexAdd()
    {}
    void generateProposalConfiguration(const Configuration&);
    FPType generateProposalPropability() throw() MTL_PURE_FUNCTION;
    /**
    If the new Configuration is accepted, this function returns it. Here the LU decomposed Matrices and the inverse matrices get updated.
    @param c the old configuration
    @return the Configuration that is proposed
    */
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used Hubbard U
    const FPType lowerpos;
    const FPType upperpos;
    const FPType contourpartlen;
    const unsigned int N;///< the number of sites
    Prng prng;///< the used PRNG
    MatType u2;
    MatType v1;
    typename Configuration::MatConf::RatioType ratio_storage;
    const FPType prefactor;
};

template<class Config>
typename GeneralMove<Config>::Configuration& ColdAtomVertexAdd<Config>::returnProposalConfiguration(const Configuration& c)
{
   //the correction vectors u2 and v1 were already calculated in generateProposalConfiguration
    this->proposalconfiguration.matcont.matrixOfOneAddedVertex(c.matcont, ratio_storage, u2, v1);
    return this->proposalconfiguration;
}

template <class Config, typename FPType>
class VertexDomain<Hubbard_Vertex<FPType>, ColdAtomVertexAdd<Config>, FPType>
{
public:
template<class PRNG>
static inline Hubbard_Vertex<FPType> generateRandomVertex(FPType lowerpos, FPType contourlen, unsigned int sites, PRNG& prng) throw()
{
 return Hubbard_Vertex<FPType>(prng.rndint(sites), lowerpos + prng.rndfloat(contourlen), (prng.template rndInteger<64>() < 32 ? DOWN : UP));//Note: SPINS SHOULD BE CHOSEN FROM -1,1...........
}
};

template<class Config>
void ColdAtomVertexAdd<Config>::generateProposalConfiguration(const Configuration& c)
{
    //Copy the old Configuration
    this->proposalconfiguration.lazyclone(c);
    //Attach a new Vertex to the proposalconfiguration
    VA_Vertex temp(VertexDomain<VA_Vertex, ColdAtomVertexAdd<Config>, FPType>::template generateRandomVertex<Prng>(lowerpos, contourpartlen, N, prng));
    this->proposalconfiguration.push_back(temp);//copy the vertex into the configuration
    //Determine the determinant and store it
    ratio_storage = c.matcont.template detOfOneAddedVertex<typename Config::GreensFunction>(temp, c, this->alpha, u2, v1);
    updateWeights(this->ratio, this->interpretableratio, this->proposalconfiguration.incphase, ratio_storage);
    return;
}

template<class Config>
typename Config::FPType ColdAtomVertexAdd<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the transition propability
    return prefactor * this->interpretableratio / static_cast<FPType>(this->proposalconfiguration.size());
}

/**
A class that encapsulates the steps for adding a Vertex
*/
template<class Config>
class VertexAdd : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< a typedef to access the Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< a typedef for the used PRNG
    typedef typename Config::Vertex VA_Vertex;///< the used Vertex Type
    typedef typename Config::FPType FPType;///< the used floating point type
    typedef typename Config::Matrix MatType;
    typedef typename Config::GreensFunction::FreeGreensFunctionReturnValueType GFType;///< the type of the Greensfunction return value
public:
    /**
    This is the Constructor for the move to add Vertices.
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _N the Number of sites
    @param a the used alpha object
    */
    VertexAdd(FPType _u, FPType cl, unsigned int _N, const Alpha<FPType>& a) :GeneralMove<Config>(1.0/2.0/*the cumulative ratio of this move*/, a),  u(_u), contourlen(cl), N(_N), prng(rand(), rand()), prefactor(u * contourlen * N)
    {}
    ~VertexAdd()
    {}
    void generateProposalConfiguration(const Configuration&);
    FPType generateProposalPropability() throw() MTL_PURE_FUNCTION;
    /**
    If the new configuration is accepted, this function returns it. Here the LU decomposed Matrices and the inverse matrices get updated.
    @param c the old configuration
    @return the Configuration that is proposed
    */
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used Hubbard U
    const FPType contourlen;///< the length of the contour on which interactions take place
    const unsigned int N;///< the number of sites
    Prng prng;///< the used PRNG
    MatType u2;
    MatType v1;
    typename Configuration::MatConf::RatioType ratio_storage;
    const FPType prefactor;
};

template<class Config>
typename GeneralMove<Config>::Configuration& VertexAdd<Config>::returnProposalConfiguration(const Configuration& c)
{
    this->proposalconfiguration.matcont.matrixOfOneAddedVertex(c.matcont, ratio_storage, u2, v1);
    return this->proposalconfiguration;
}

template <class Config, typename FPType>
class VertexDomain<Hubbard_Vertex<FPType>, VertexAdd<Config>, FPType>
{
public:
template<class PRNG>
static inline Hubbard_Vertex<FPType> generateRandomVertex(FPType contourlen, unsigned int sites, PRNG& prng) throw()
{
return Hubbard_Vertex<FPType>(prng.rndint(sites), prng.rndfloat(contourlen),
                      (prng.template rndInteger<64>() < 32 ? DOWN : UP));//Spins should be chosen from {-1, 1}, to be in accordance with Fakhers paper
}
};

template<class Config>
void VertexAdd<Config>::generateProposalConfiguration(const Configuration& c)
{
    //Copy the old Configuration
    this->proposalconfiguration.lazyclone(c);
    //Attach a new Vertex to the proposalconfiguration
    VA_Vertex temp(VertexDomain<VA_Vertex, VertexAdd<Config>, FPType>::template generateRandomVertex<Prng>(contourlen, N, prng));
    this->proposalconfiguration.push_back(temp);//copy the vertex into the configuration
    //Determine the determinant and store it
    ratio_storage = c.matcont.template detOfOneAddedVertex<typename Config::GreensFunction>(temp, c, this->alpha, u2, v1);
    updateWeights(this->ratio, this->interpretableratio, this->proposalconfiguration.incphase, ratio_storage);
    return;
}

template<class Config>
typename Config::FPType VertexAdd<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the transition propability
    return prefactor * this->interpretableratio / static_cast<FPType>(this->proposalconfiguration.size());
}


/**
A class that encapsulates the steps for adding a Vertex
*/
template<class Config>
class MTMVertexAdd : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< a typedef to access the Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< a typedef for the used PRNG
    typedef typename Config::Vertex VA_Vertex;///< the used Vertex Type
    typedef typename Config::FPType FPType;///< the used floating point type
    typedef typename Config::Matrix MatType;
    typedef typename Config::GreensFunction::FreeGreensFunctionReturnValueType GFType;///< the type of the Greensfunction return value
public:
    /**
    This is the Constructor for the move to add Vertices.
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _N the Number of sites
    @param a the used alpha object
    */
    MTMVertexAdd(FPType _u, FPType cl, unsigned int _N, const Alpha<FPType>& a) :GeneralMove<Config>(1.0/2.0/*the cumulative ratio of this move*/, a),  u(_u), contourlen(cl), N(_N), prng(rand(), rand()), prefactor(u * contourlen * N)
    {
    }
    ~MTMVertexAdd()
    {
    }
    void generateProposalConfiguration(const Configuration&);
    FPType generateProposalPropability() throw() MTL_PURE_FUNCTION;
    /**
    If the new configuration is accepted, this function returns it. Here the LU decomposed Matrices and the inverse matrices get updated.
    @param c the old configuration
    @return the Configuration that is proposed
    */
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used Hubbard U
    const FPType contourlen;///< the length of the contour on which interactions take place
    const unsigned int N;///< the number of sites
    Prng prng;///< the used PRNG
    MatType u2;
    MatType v1;
    typename Configuration::MatConf::RatioType ratio_storage;
    const FPType prefactor;
    FPType a1;//the sum of all added vertices
    FPType b1;//the sum of all removed vertices
};

template<class Config>
typename GeneralMove<Config>::Configuration& MTMVertexAdd<Config>::returnProposalConfiguration(const Configuration& c)
{
    this->proposalconfiguration.matcont.matrixOfOneAddedVertex(c.matcont, ratio_storage, u2, v1);
    return this->proposalconfiguration;
}

template <class Config, typename FPType>
class VertexDomain<Hubbard_Vertex<FPType>, MTMVertexAdd<Config>, FPType>
{
public:
template<class PRNG>
static inline Hubbard_Vertex<FPType> generateRandomVertex(FPType contourlen, unsigned int sites, PRNG& prng) throw()
{
return Hubbard_Vertex<FPType>(prng.rndint(sites), prng.rndfloat(contourlen),
                      (prng.template rndInteger<64>() < 32 ? DOWN : UP));//Spins should be chosen from {-1, 1}, to be in accordance with Fakhers paper
}
};

template<class Config>
void MTMVertexAdd<Config>::generateProposalConfiguration(const Configuration& c)
{
  //generate all the proposal Vertices
  std::vector<VA_Vertex> nv;
  std::cout<<"old order: "<<c.size()<<std::endl;
  for(uint k = 0; k < c.size(); ++k)
    nv.push_back(VertexDomain<VA_Vertex, VertexAdd<Config>, FPType>::template generateRandomVertex<Prng>(contourlen, N, prng));
  //generate the reference states.
  uint* indices = new uint[c.size()+1];
  memset(indices, 0, sizeof(uint) * (c.size() + 1));
  for(uint k = 0; k < (c.size()-1); ++k)
    indices[prng.rndint(c.size() + 1)]++;
  //determine which Vertex to choose. We will also call the final Vertex y
  FPType* rs = new FPType[nv.size()];
  for(uint k = 0; k < nv.size(); ++k)
  {
    auto temp = c.matcont.template detOfOneAddedVertex<typename Config::GreensFunction>(nv[k], c, this->alpha, u2, v1);
    rs[k] = -temp.first * temp.second;
//    std::cout<<rs[k]<<std::endl;
  }
  FPType nomin = 0;//The norm will contain the sum of all proposed vertices that we need later.
  for(uint k = 0; k < nv.size(); ++k)
    nomin += rs[k];
  std::cout<<nomin<<std::endl;
  for(uint k = 0; k < nv.size(); ++k)
  {
    rs[k] /= nomin;
  }
  FPType avw = 1.0/c.size();
//   for(uint k = 0; k < nv.size(); ++k)
//   {
//    std::cout<<k<<" "<<-std::log(rs[k])<<" "<<nv[k]<<std::endl;
//   }
  
  FPType rndnr = prng.rndfloat();
  //check which index we will choose
  uint yidx = 0;
  FPType inc = 0;
  for(; yidx < nv.size() && inc < rndnr; ++yidx)
  {
    inc += rs[yidx];
//    std::cout<<yidx<<" "<<inc<<std::endl;
  }
  yidx = yidx - 1;
  //now yidx contains the index of the vertex that we propose to move to.
 std::cout<<rndnr<<" -> "<<yidx<<std::endl;
 //Copy the old Configuration
 this->proposalconfiguration.lazyclone(c);
 this->proposalconfiguration.push_back(nv[yidx]);
 //Determine the determinant and store it
 ratio_storage = c.matcont.template detOfOneAddedVertex<typename Config::GreensFunction>(nv[yidx], c, this->alpha, u2, v1);
 this->proposalconfiguration.matcont.matrixOfOneAddedVertex(c.matcont, ratio_storage, u2, v1);
 //Now we determine the total weight
 b1 = 1.0;//The contribution of going to the old contribution which has to be present
 for(uint k = 0; k < c.size() + 1; ++k)
 {
   if(indices[k] > 0)
   {
     auto temp = this->proposalconfiguration.matcont.detOfOneRemovedVertex(k);
     b1 += indices[k] * (-temp.first * temp.second);
   }
 }
 updateWeights(this->ratio, this->interpretableratio, this->proposalconfiguration.incphase, ratio_storage);
 a1 = nomin;
    return;
}

template<class Config>
typename Config::FPType MTMVertexAdd<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the transition propability
    return prefactor * this->interpretableratio *a1/ static_cast<FPType>(this->proposalconfiguration.size() * b1);
}

/**
A class that encapsulates the steps for moving a Vertex in time
*/
template<class Config>
class MoveVertexinTime : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< a typedef to access the Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< a typedef for the used PRNG
    typedef typename Config::Vertex VA_Vertex;///< the used Vertex Type
    typedef typename Config::FPType FPType;///< the used floating point type
public:
    /**
    This is the Constructor for the move to add Vertices.
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _N the Number of sites
    @param a the used alpha object
    */
    MoveVertexinTime(float _u, float _beta, unsigned int _N, const Alpha<FPType>& a) :GeneralMove<Config>(3.0/7.0/*the cumulative ratio of this move*/, a),  u(_u), beta(_beta), N(_N), prng(rand(), static_cast<float>(rand()))
    {}
    ~MoveVertexinTime()
    {}
    void generateProposalConfiguration(const Configuration&);
    typename Config::FPType generateProposalPropability() throw();
    /**
    If the new Configuration is accepted, this function returns it. Here the LU decomposed Matrices and the inverse matrices get updated.
    @param c the old configuration
    @return the Configuration that is proposed
    */
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used Hubbard U
    const FPType beta;///< the used inverse temperature
    const unsigned int N;///< the number of sites
    Prng prng;///< the used PRNG
};

template<class Config>
typename GeneralMove<Config>::Configuration& MoveVertexinTime<Config>::returnProposalConfiguration(const Configuration& c)
{
    return this->proposalconfiguration;
}

template<class Config>
void MoveVertexinTime<Config>::generateProposalConfiguration(const Configuration& c)
{
    this->oldcfg = c;
    this->proposalconfiguration = c;
    this->proposalconfiguration.lazyreset();
    for (typename Configuration::iterator it = this->proposalconfiguration.begin(); it != this->proposalconfiguration.end(); ++it)
        it->tau = prng.rndfloat(beta);
    this->proposalconfiguration.template reevaluate<typename Config::GreensFunction>(this->alpha);
    return;
}

template<class Config>
typename Config::FPType MoveVertexinTime<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the transition propability
    FPType denominator = this->oldcfg.interpretableweight;

    //Now do the same for the nominator
    //Calculate the Nominator of the transition propability
    FPType nominator = this->proposalconfiguration.interpretableweight;
    return nominator / denominator;//Calculate ratio
}


/**
A class that encapsulates the steps for moving the vertex
*/
template<class Config>
class MoveVertexinSpace : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< a typedef to access the Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< a typedef for the used PRNG
    typedef typename Config::Vertex VA_Vertex;///< the used Vertex Type
    typedef typename Config::FPType FPType;///< the used floating point type
public:
    /**
    This is the Constructor for the move to add Vertices.
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _N the Number of sites
    @param a the used alpha object
    */
    MoveVertexinSpace(float _u, float _beta, unsigned int _N, const Alpha<FPType>& a) :GeneralMove<Config>(4.0/7.0/*the cumulative ratio of this move*/, a),  u(_u), beta(_beta), N(_N), prng(/*rand()*/160487425, static_cast<float>(/*rand()*/ 27957287))
    {}
    ~MoveVertexinSpace()
    {}
    void generateProposalConfiguration(const Configuration&);
    typename Config::FPType generateProposalPropability() throw();
    /**
    If the new Configuration is accepted, this function returns it. Here the LU decomposed Matrices and the inverse matrices get updated.
    @param c the old configuration
    @return the Configuration that is proposed
    */
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used Hubbard U
    const FPType beta;///< the used inverse temperature
    const unsigned int N;///< the number of sites
    Prng prng;///< the used PRNG
};

template<class Config>
typename GeneralMove<Config>::Configuration& MoveVertexinSpace<Config>::returnProposalConfiguration(const Configuration& c)
{
    return this->proposalconfiguration;
}

template<class Config>
void MoveVertexinSpace<Config>::generateProposalConfiguration(const Configuration& c)
{
    this->oldcfg = c;
    this->proposalconfiguration = c;
    this->proposalconfiguration.lazyreset();
    for (typename Configuration::iterator it = this->proposalconfiguration.begin(); it != this->proposalconfiguration.end(); ++it)
        it->site = prng.rndint(2);
    this->proposalconfiguration.template reevaluate<typename Config::GreensFunction>(this->alpha);
    return;
}

template<class Config>
typename Config::FPType MoveVertexinSpace<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the transition propability
    FPType denominator = this->oldcfg.interpretableweight;

    //Now do the same for the nominator
    //Calculate the Nominator of the transition propability
    FPType nominator = this->proposalconfiguration.interpretableweight;
    return nominator / denominator;//Calculate ratio
}

/**
Some Pseudo - Move, which could mean a spin-flip
*/
template<class Config>
class SpinFlip : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;
    typedef typename GeneralMove<Config>::Prng Prng;
    typedef typename Config::FPType FPType;///< the used floating point type
public:
    /**
    This is the constructor for the move to flip the spins in the configuration
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _n the Number of sites
    @param a the used alpha Object
    */
    SpinFlip(FPType _u, FPType _beta, const Alpha<FPType>& a) : GeneralMove<Config>(1.0, a), u(_u), prng(rand(), rand())
    {
    }
    ~SpinFlip() {}
    void generateProposalConfiguration(const Configuration&);
    FPType generateProposalPropability() throw()
    {
      return this->interpretableratio;
    }
    Configuration& returnProposalConfiguration(const Configuration& c)
    {
      if(c.size() > 0)
      {
	this->proposalconfiguration.matcont.matrix_of_a_flipped_Ising_Spin(flip_pos, c.matcont, ratio_storage, gstore);
      }
      return this->proposalconfiguration;
    }
private:
    FPType u;
    Prng prng;
    uint flip_pos;
    typename Configuration::MatConf::RatioType ratio_storage;
    typename Configuration::MatConf::RatioType gstore;
};

template<class Config>
void SpinFlip<Config>::generateProposalConfiguration(const Configuration& c)
{
    if(c.size() > 0)
    {
      this->proposalconfiguration.lazyclone(c);
      flip_pos = prng.rndint(c.size());
      this->proposalconfiguration[flip_pos].spin = !c[flip_pos].spin;
      ratio_storage = c.matcont.template det_of_a_flippedIsingSpin<Configuration>(flip_pos, c, this->alpha, gstore);
      updateWeights(this->ratio, this->interpretableratio, this->proposalconfiguration.incphase, ratio_storage);
    }
    else
    {
      this->proposalconfiguration = c;
      this->interpretableratio = 1.0;
    }
    return;
}

template<class Config>
class AddTwoVertices : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< a typedef to access the Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< a typedef for the used PRNG
    typedef typename Config::Vertex VA_Vertex;///< the used Vertex Type
    typedef typename Config::FPType FPType;///< the used floating point type
public:
    /**
    This is the Constructor for the move to add Vertices.
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _N the Number of sites
    @param a the used alpha object
    */
    AddTwoVertices(float _u, float _beta, unsigned int _N, const Alpha<FPType>& a) :GeneralMove<Config>(4.0/5.0/*the cumulative ratio of this move*/, a),  u(_u), beta(_beta), N(_N), prng(/*rand()*/160487425, static_cast<float>(/*rand()*/ 27957287))
    {}
    ~AddTwoVertices()
    {}
    void generateProposalConfiguration(const Configuration&);
    typename Config::FPType generateProposalPropability() throw();
    /**
    If the new Configuration is accepted, this function returns it. Here the LU decomposed Matrices and the inverse matrices get updated.
    @param c the old configuration
    @return the Configuration that is proposed
    */
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used Hubbard U
    const FPType beta;///< the used inverse temperature
    const unsigned int N;///< the number of sites
    Prng prng;///< the used PRNG
};

template<class Config>
typename GeneralMove<Config>::Configuration& AddTwoVertices<Config>::returnProposalConfiguration(const Configuration& c)
{
    return this->proposalconfiguration;
}

template<class Config>
void AddTwoVertices<Config>::generateProposalConfiguration(const Configuration& c)
{
    this->oldcfg = c;
    this->proposalconfiguration = c;
    this->proposalconfiguration.lazyreset();
    int s = prng.template rndInteger<32>();//to base the decision of the Spin not only on the least significant bit
    s = (s < 16 ? 1 : -1);//Spins should be chosen from {-1, 1}, to be in accordance with Fakhers paper
    unsigned int n = (prng.rndfloat() < 0.5 ? 0:1);
    this->proposalconfiguration.push_back( VA_Vertex(prng.rndfloat(beta), n, s));
    s = prng.template rndInteger<32>();//to base the decision of the Spin not only on the least significant bit
    s = (s < 16 ? 1 : -1);//Spins should be chosen from {-1, 1}, to be in accordance with Fakhers paper
    n = (prng.rndfloat() < 0.5 ? 0:1);
    this->proposalconfiguration.push_back( VA_Vertex(prng.rndfloat(beta), n, s));
    this->proposalconfiguration.template reevaluate<typename Config::GreensFunction>(this->alpha);
    return;
}

template<class Config>
typename Config::FPType AddTwoVertices<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the transition propability
    FPType denominator = static_cast<FPType>(this->proposalconfiguration.size() * (this->proposalconfiguration.size()-1)) * this->oldcfg.interpretableweight;

    //Now do the same for the nominator
    //Calculate the Nominator of the transition propability
    FPType nominator = u*u * beta*beta * N*N * this->proposalconfiguration.interpretableweight;
    return nominator / denominator;//Calculate ratio
}


template<class Config>
class RemoveTwoVertices : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< a typedef to access the Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< a typedef for the used PRNG
    typedef typename Config::Vertex VA_Vertex;///< the used Vertex Type
    typedef typename Config::FPType FPType;///< the used floating point type
public:
    /**
    This is the Constructor for the move to add Vertices.
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _N the Number of sites
    @param a the used alpha object
    */
    RemoveTwoVertices(float _u, float _beta, unsigned int _N, const Alpha<FPType>& a) :GeneralMove<Config>(5.0/5.0/*the cumulative ratio of this move*/, a),  u(_u), beta(_beta), N(_N), prng(rand(), static_cast<float>(rand()))
    {}
    ~RemoveTwoVertices()
    {}
    void generateProposalConfiguration(const Configuration&);
    typename Config::FPType generateProposalPropability() throw();
    /**
    If the new Configuration is accepted, this function returns it. Here the LU decomposed Matrices and the inverse matrices get updated.
    @param c the old configuration
    @return the Configuration that is proposed
    */
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used Hubbard U
    const FPType beta;///< the used inverse temperature
    const unsigned int N;///< the number of sites
    Prng prng;///< the used PRNG
};

template<class Config>
typename GeneralMove<Config>::Configuration& RemoveTwoVertices<Config>::returnProposalConfiguration(const Configuration& c)
{
    return this->proposalconfiguration;
}

template<class Config>
void RemoveTwoVertices<Config>::generateProposalConfiguration(const Configuration& c)
{
    this->oldcfg = c;
    this->proposalconfiguration = c;
    this->proposalconfiguration.lazyreset();
    const unsigned int order = c.size();

    if (order < 3)
    {
        this->proposalconfiguration.resetmatrices();
        if (order == 0) return;//Nothing to do without vertex.
        else
        {
            //if we are at the first order the next configuration will have zero vertices and thus have determinants and a weight of 1.0
            this->proposalconfiguration.pop_back();
            this->proposalconfiguration.pop_back();
            return;
        }
    }
    //Generate randomly the vertex position which we remove
    unsigned int tmppos = prng.rndint(order);
    //Erase the Vertex
    typename Configuration::iterator it = this->proposalconfiguration.begin();
    advance(it, tmppos);
    this->proposalconfiguration.erase(it);
    tmppos = prng.rndint(order);
    //Erase the Vertex
    it = this->proposalconfiguration.begin();
    advance(it, tmppos);
    this->proposalconfiguration.erase(it);
    this->proposalconfiguration.template reevaluate<typename Config::GreensFunction>(this->alpha);
    return;
}

template<class Config>
typename Config::FPType RemoveTwoVertices<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the Transition Propability
    FPType nominator = static_cast<FPType>(this->oldcfg.size() * (this->oldcfg.size()-1) ) * this->proposalconfiguration.interpretableweight;
    //Now do the same for the nominator
    //Calculate the Nominator of the Transition Propability
    FPType denominator = u*u * beta*beta * N*N * this->oldcfg.interpretableweight;
    return nominator / denominator;//Calculate ratio
}

/**
A class for removing a Vertex
*/
template<class Config>
class VertexRemove : public GeneralMove<Config>
{
    typedef typename GeneralMove<Config>::Configuration Configuration;///< the type of the used Configuration
    typedef typename GeneralMove<Config>::Prng Prng;///< the Type of the used PRNG
    typedef typename Config::Vertex VR_Vertex;///< the Type of the used Vertex
    typedef typename Config::FPType FPType;///< the floating point data type to use for the calculations
    typedef typename Config::GreensFunction::FreeGreensFunctionReturnValueType GFType;
    typedef typename Config::Matrix Matrix;
public:
    /**
    This is the Constructor for the move to remove Vertices.
    @param _u the Hubbard U
    @param _beta the inverse temperature
    @param _N the Number of sites
    @param a the used alpha object
    */
    VertexRemove(FPType _u, FPType _beta, unsigned int _n, const Alpha<FPType>& a) : GeneralMove<Config>(1.0/*the cumulative ratio of this move*/, a), u(_u), beta(_beta), N(_n), prefactor(1.0/(u * beta * static_cast<FPType>(N))), prng(rand(),rand())
    {}
    void generateProposalConfiguration(const Configuration&);
    FPType generateProposalPropability() throw() MTL_PURE_FUNCTION;
    Configuration& returnProposalConfiguration(const Configuration& c);
private:
    const FPType u;///< the used parameter of the Hubbard Interaction
    const FPType beta;///< the used inverse temperature beta
    const unsigned int N;///< the number of sites
    const FPType prefactor;
    Prng prng;///< the used PRNG
    unsigned int tmppos;//internal state don't mess around with it....
};

template<class Config>
typename Config::FPType VertexRemove<Config>::generateProposalPropability() throw()
{
    //Calculate the denominator of the Transition Propability
//    std::cout<<"[VertexRemove] this->interpretableratio = "<<this->interpretableratio<<std::endl; 
//    std::cout<<"[VertexRemove] this->proposalconfiguration.size() = "<<this->proposalconfiguration.size()<<std::endl;
//    std::cout<<"[VertexRemove] prefactor = "<<prefactor<<std::endl;
    return this->interpretableratio * static_cast<FPType>(this->proposalconfiguration.size() + 1) * prefactor;
}

template <class Config>
void VertexRemove<Config>::generateProposalConfiguration(const Configuration& c)
{
    //Copy the old configuration
    this->proposalconfiguration.lazyclone(c);
    //reset all weight related things to initial values
    const unsigned int order = static_cast<unsigned int>(c.size());
    if (order < 2)
    {
        this->proposalconfiguration.resetmatrices();
        if (order == 0)
        {
            this->proposalconfiguration.incphase = this->ratio = this->interpretableratio = 1.0;
        }
        if (order == 1)
        {
            this->proposalconfiguration.pop_back();
            updateWeights(this->ratio, this->interpretableratio, this->proposalconfiguration.incphase, c.matcont.detOfOneRemovedVertex(0));
            this->proposalconfiguration.incphase = 1.0;
        }
        return;
    }
    //Generate randomly the vertex position which we remove
    tmppos = prng.rndint(order);
    //Erase the Vertex
    typename Configuration::iterator it = this->proposalconfiguration.begin();
    advance(it, tmppos);
    //A new Idea: Since the order of the vertices doesn't matter we swap the to-be-removed element to the end and then invalidate the last position
    iter_swap(it, this->proposalconfiguration.end()-1);
    this->proposalconfiguration.pop_back();
    //Note that the new determinant(the one with one vertex less) are related by(A similar relation holds for blockwise updates):
//    det(A_{n-1}) = det(A_{n}) * (A^{-1})_{tmppos, tmppos}
//    std::cout<<this->proposalconfiguration<<std::endl;
    updateWeights(this->ratio, this->interpretableratio, this->proposalconfiguration.incphase, c.matcont.detOfOneRemovedVertex(tmppos) );
    return;
}

template<class Config>
typename GeneralMove<Config>::Configuration& VertexRemove<Config>::returnProposalConfiguration(const Configuration& conf)
{
    this->proposalconfiguration.matcont.template matrixOfOneRemovedVertex<typename Config::GreensFunction> (conf, this->tmppos, this->alpha, this->proposalconfiguration);
//    std::cout<<"[Remove]: Inverse of inverse... :"<<MTLICCL::inverse(this->proposalconfiguration.matcont.mat)<<std::endl;
    return this->proposalconfiguration;
}

/**
Structure for the DDQMC Modell
*/
template <class _Config, class Moves>
class DDQMC
{
public:
    typedef typename _Config::Container Configuration;///< the type of the used container
    typedef _Config Config;///< a typedef for the Configuration
    typedef typename Config::FPType FPType;///< the type used for floating point calculations
    GeneralMove<_Config>* moves[TListLength<Moves>::Ret];///< The array with all supported moves
    /** Allocates the moves
    @param alpha the alpha object that the moves should use
    */
    inline DDQMC(const Alpha<FPType>& alpha, const Parameters&);
    /**
    Destructor which frees the memory used by the model and by the Greensfunction tables
    */
    inline ~DDQMC();
private:
};

template <class _Config, class Moves>
DDQMC<_Config, Moves>::~DDQMC()
{
    //free the memory
    for(int k = 0; k < TListLength<Moves>::Ret; ++k)
        delete moves[k];
    // free the memory used by the Greensfunction tables
    return;
}

template <class Greensfunction, class Moves>
struct DDQMC_Init_Helper
{
template <class BasicMove>
static inline void init(BasicMove *const, const Alpha<typename BasicMove::_Config::FPType>& alpha, const Parameters& curparams);
};

template <class _Config, class Moves>
DDQMC<_Config, Moves>::DDQMC(const Alpha<FPType>& alpha, const Parameters& curparams)
{
    //Get the model parameters from the config files
    DDQMC_Init_Helper<typename _Config::GreensFunction, Moves>::init(moves, alpha, curparams);
    //Some Selftest for the SIAM Greensfunctiom
    typedef typename _Config::Vertex Vertex;
/*    Vertex stuff1[25] ={
    Vertex(36.203114, -1),
    Vertex( 0.571121,  1),
    Vertex(38.499023, -1),
    Vertex(21.000648,  1),
    Vertex( 2.295501,  1),
    Vertex(24.356186,  1),
    Vertex(16.662104, -1),
    Vertex( 2.053399,  1),
    Vertex(29.232674, -1),
    Vertex( 3.354945,  1),
    Vertex(31.139961, -1),
    Vertex(38.993202,  1),
    Vertex(24.845573,  1),
    Vertex(11.577677,  1),
    Vertex(35.062347, -1),
    Vertex(13.083897,  1),
    Vertex(38.571632,  1),
    Vertex(39.388660,  1),
    Vertex(35.063751, -1),
    Vertex(30.486607, -1),
    Vertex(19.694229,  1),
    Vertex(39.263580,  1),
    Vertex(23.488436,  1),
    Vertex(39.819214,  1),
    Vertex(16.046282, -1)
    };
    Vertex stuff2[9] = 
    {
    Vertex(16.740049,  1),
    Vertex(27.486649,  1),
    Vertex(22.019007,  1),
    Vertex( 6.477389, -1),
    Vertex(32.903568, -1),
    Vertex( 0.559886, -1),
    Vertex( 5.072774,  1),
    Vertex( 7.599339,  1),
    Vertex(10.816016, -1)
    };
    Vertex stuff3[2] =
    {
    Vertex(13.184099, -1),
    Vertex(27.140181,  1)
    };
    Vertex stuff4[28] = {
    Vertex(35.497265, -1), Vertex(36.457531, -1), Vertex(31.823751, 1), Vertex(0.571121, 1), Vertex(15.235376, 1), Vertex(3.165603, 1), Vertex(13.516393, 1), Vertex(24.356186, 1), Vertex(27.856997, -1), Vertex(16.662104, -1), Vertex(1.496163, 1), Vertex(29.232674, -1), Vertex(3.354945, 1), Vertex(31.139961, -1), Vertex(38.993202, 1), Vertex(33.855698, -1), Vertex(35.062347, -1), Vertex(19.885864, -1), Vertex(13.083897, 1), Vertex(38.571632, 1), Vertex(39.388660, 1), Vertex(8.378150, -1), Vertex(35.063751, -1), Vertex(30.486607, -1), Vertex(39.263580, 1), Vertex(23.488436, 1), Vertex(39.819214, 1), Vertex(16.046282, -1)
    };
    typename _Config::Container cont;
    for(int k = 0; k < 28; ++k)
    cont.push_back(stuff4[k]);
    Reevaluate<typename _Config::DeterminantType>::template reevaluate<
    typename _Config::GreensFunction, Configuration>(cont,alpha);
  */  
/*    std::ofstream file_re("GFTestdownup.dat");
    std::ofstream file_im("GFTestdowndown_im.dat");
    Vertex v1(0, 0.0, DOWN);
    for(float tau = 0.0; tau <= curparams.beta; tau += 0.001)
    {
    Vertex v2(1, tau, UP);
    typename _Config::GreensFunction::FreeGreensFunctionReturnValueType temp = _Config::GreensFunction::eval(v2, v1);
    file_re<<tau<<" "<<temp<<std::endl;
    file_im<<tau<<" "<<imag(temp)<<std::endl;
    }
    exit(-1);*/
/*    std::cout.precision(10);
    std::cout<<cont<<std::endl;
    std::cout<<cont.up.inverse<<endl;
    std::cout<<cont.down.inverse<<endl;
    for(int k = 0; k < cont.up.inverse.Rows(); ++k)
    for(int i = 0; i < cont.up.inverse.Rows(); ++i)
    std::cout<<cont.up.inverse(k,i)<<std::endl;;
    exit(-1);*/
    /*
    some selftest for debugging the Realtime Hubbard Greensfunction
    typedef typename _Config::Vertex Vertex;
    Vertex v1(0.1, 0, UP);
    Vertex v2(0.2, 1, UP);
    Vertex v3(0.9, 0, UP);
    Vertex v4(0.8, 1, UP);
    Vertex v5(5.0, 1, UP);
    std::cout<<_Config::GreensFunction::eval(v1, v1)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v2, v2)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v3, v3)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v4, v4)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v5, v5)<<std::endl;

    std::cout<<_Config::GreensFunction::eval(v1, v2)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v1, v3)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v1, v4)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v1, v5)<<std::endl;

    std::cout<<_Config::GreensFunction::eval(v2, v1)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v2, v3)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v2, v4)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v2, v5)<<std::endl;

    std::cout<<_Config::GreensFunction::eval(v3, v1)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v3, v2)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v3, v4)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v3, v5)<<std::endl;

    std::cout<<_Config::GreensFunction::eval(v4, v1)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v4, v2)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v4, v3)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v4, v5)<<std::endl;

    std::cout<<_Config::GreensFunction::eval(v5, v1)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v5, v2)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v5, v3)<<std::endl;
    std::cout<<_Config::GreensFunction::eval(v5, v4)<<std::endl;
    exit(-1);*/
    return;
}
#endif
