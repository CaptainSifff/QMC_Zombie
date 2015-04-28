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
#ifndef OBSERVABLES_H
#define OBSERVABLES_H
#include "AverageSign.h"
#include "ObservableBase.h"
#include "ObservableContainer.h"
#include "Parameters.h"
#include "Greensfunction.h"
#include "libFourier.h"
#include <valarray>

#ifdef _OPENMP
#include <omp.h>
#endif

#define TWOPI 2.0*M_PIl

template <class C, class S, class Communication>
struct Observable_Parameters
{
    typedef C Configuration;
    typedef S SignType;
    typedef Communication Comm;
};

template <class Configuration>
void densitydensityCorrelation_dry(DryRun<typename Configuration::value_type, typename Configuration::DetType>& func,
                                   int site_i, const typename DryRun<typename Configuration::value_type, typename Configuration::DetType>::FPType s_i, int spin_i,
                                   int site_j, const typename DryRun<typename Configuration::value_type, typename Configuration::DetType>::FPType s_j, int spin_j)
{
    const typename DryRun<typename Configuration::value_type, typename Configuration::DetType>::FPType tiny = 0.00000001;
    if (spin_i == UP)
    {
        if (spin_j == DOWN)
        {
            //UP-DOWN
            func.template onSector<UP>(site_i, s_i, site_i, s_i);
            func.template onSector<DOWN>(site_j, s_j, site_j, s_j);
        }
        else
        {
            //UP-UP
            if ((site_i != site_j) || !fpequal(s_i, s_j))
            {
                func.template onSector<UP>(site_i, s_i, site_i, s_i);
                func.template onSector<UP>(site_j, s_j, site_j, s_j);
                func.template onSector<UP>(site_j, s_j, site_i, s_i);
                func.template onSector<UP>(site_i, s_i, site_j, s_j);
            }
            else
            {
                func.template onSector<UP>(site_i, s_i, site_i, s_i);
            }
        }
    }
    else
    {
        if (spin_j == DOWN)
        {
            if ((site_i != site_j) || !fpequal(s_i, s_j))
            {
                //DOWN-DOWN
                func.template onSector<DOWN>(site_i, s_i, site_i, s_i);
                func.template onSector<DOWN>(site_j, s_j, site_j, s_j);
                func.template onSector<DOWN>(site_j, s_j, site_i, s_i);
                func.template onSector<DOWN>(site_i, s_i, site_j, s_j);
            }
            else
            {
                func.template onSector<DOWN>(site_i, s_i, site_i, s_i);
            }
        }
        else
        {
            //DOWN - UP
            func.template onSector<DOWN>(site_i, s_i, site_i, s_i);
            func.template onSector<UP>(site_j, s_j, site_j, s_j);
        }
    }
    return;
}

template <class Configuration>
inline typename Configuration::DetType densitydensityCorrelation(const DoWick<typename Configuration::value_type, typename Configuration::DetType>& dowick, int site_i, typename DoWick<typename Configuration::value_type, typename Configuration::DetType>::FPType s_i, int spin_i, int site_j, typename DoWick<typename Configuration::value_type, typename Configuration::DetType>::FPType s_j, int spin_j)
{
    typename Configuration::DetType retval;
    const typename DoWick<typename Configuration::value_type, typename Configuration::DetType>::FPType tiny = 0.00000001;
    if (spin_i == UP)
    {
        if (spin_j == DOWN)
        {
            //UP-DOWN
            retval = dowick.template onSector<UP>(site_i, s_i, site_i, s_i)
                     * dowick.template onSector<DOWN>(site_j, s_j, site_j, s_j);
        }
        else
        {
            //UP-UP
            if ((site_i != site_j) || !fpequal(s_i, s_j))
            {
                retval = (dowick.template onSector<UP>(site_i, s_i, site_i, s_i)
                          * dowick.template onSector<UP>(site_j, s_j, site_j, s_j)
                          - dowick.template onSector<UP>(site_j, s_j, site_i, s_i)
                          * dowick.template onSector<UP>(site_i, s_i, site_j, s_j));
            }
            else
            {//special care for equaltime greensfunction
                retval = dowick.template onSector<UP>(site_i, s_i, site_i, s_i);
            }
        }
    }
    else
    {
        if (spin_j == DOWN)
        {
            if ((site_i != site_j) || !fpequal(s_i, s_j))
            {
                //DOWN-DOWN
                retval = (dowick.template onSector<DOWN>(site_i, s_i, site_i, s_i)
                          * dowick.template onSector<DOWN>(site_j, s_j, site_j, s_j)
                          - dowick.template onSector<DOWN>(site_j, s_j, site_i, s_i)
                          * dowick.template onSector<DOWN>(site_i, s_i, site_j, s_j)
                         );
            }
            else
            {//special care for equaltime greensfunction
                retval = dowick.template onSector<DOWN>(site_i, s_i, site_i, s_i);
            }
        }
        else
        {
            //DOWN - UP
            retval = dowick.template onSector<DOWN>(site_i, s_i, site_i, s_i)
                     * dowick.template onSector<UP>(site_j, s_j, site_j, s_j);
        }
    }
    return retval;
}

/**
A "generic" twoparticle Greensfunction of this form:
G = <c^dagger_a c_b c^dagger_c c_d >
*/
template <class Configuration>
void genericTwoParticleGreensfunction_dry(DryRun<typename Configuration::value_type, typename Configuration::DetType>& func,
const typename Configuration::value_type& va,
const typename Configuration::value_type& vb,
const typename Configuration::value_type& vc,
const typename Configuration::value_type& vd)
{
  if((vb == vc) && (va == vb) && (vc == vd))
    func(va, va);
  else
  {
    func(va, vb);
    func(va, vd);
    func(vc, vb);
    func(vc, vd);
  }
    return;
}

/**
A "generic" twoparticle Greensfunction:
It is like that:
G = <c^dagger_a c_b c^dagger_c c_d >
*/
template <class Configuration>
typename Configuration::DetType genericTwoParticleGreensfunction(const DoWick<typename Configuration::value_type, typename Configuration::DetType>& dowick,
typename Configuration::value_type va,
typename Configuration::value_type vb,
typename Configuration::value_type vc,
typename Configuration::value_type vd)
{
  typename Configuration::DetType retval;
  if((vb == vc) && (va == vb) && (vc == vd))
  {
    retval = dowick(va, va);
  }
  else
  {
    retval =  dowick(va, vb) * dowick(vc, vd) - dowick(va, vd) * dowick(vc, vb);
  }
  return retval;
}

/**
A class for measuring the average Order.
*/
template <class Config>
class AverageOrder : public Network_Cache<Config, typename Config::SignType>
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Config::Configuration::FPType FPType;
    typedef typename Config::SignType SignType;
    typedef typename Config::SignType ObservableType;///< the AverageOrder is essentially a floating Point type, but during calculations it can be complex
    /**
    The Constructor for the Average Order
    */
    AverageOrder(typename Config::Comm& n) throw() : Network_Cache<Config, typename Config::SignType>(n, "AverageOrder"), averageorder(0.0), configurationLength(0.0)
    {
    }
    inline void dryrun(DryRun<typename Configuration::value_type, SignType>&) {}
    /**
    This determines the Average Order for the given configuration
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, SignType>&);
private:
    SignType averageorder;///< here we store the physical AverageOrder
    FPType configurationLength;///< this is the real Average Order of the Data Structure
};

template <class Config>
void AverageOrder<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, SignType>&)
{
    FPType currentorder = static_cast<FPType>(configuration.size());
    //now update the average order
    configurationLength += currentorder;
    SignType obs = currentorder * configuration.phase;//as for every other observable, account for the sign of the Configuration
    averageorder += obs;
    this->add_bin(obs);
    return;
}

/**
A class for measuring the ParticleNumber at a certain hardcoded time
*/
template <class Config>
class ParticleNumber : public Network_Cache<Config, typename Config::SignType>
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType SignType;
    typedef SignType ObservableType;///< the ParticleNumber is essentially a floating point type, but during calculations it can be complex
    /**
    The Constructor for the Average Order
    */
    ParticleNumber(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, typename Config::SignType>(n, "ParticleNumber"), t_M(params.t_exp), len(params.N)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, SignType>&);
    /**
    This determines the Average Order for the given Configuration
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, SignType>&);
private:
    const double& t_M;
    const uint32_t& len;
};

template <class Config>
void ParticleNumber<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, SignType>& func)
{
    ObservableType obs(func.template onSector<UP>(0, t_M, 0, t_M));
    obs += func.template onSector<DOWN>(0, t_M, 0, t_M);
    for (unsigned int k = 1; k < len; ++k)
    {
        obs += func.template onSector<UP>(k, t_M, k, t_M);
        obs += func.template onSector<DOWN>(k, t_M, k, t_M);
    }
    obs *= configuration.phase;
    this->add_bin(obs);
    return;
}

template <class Config>
void ParticleNumber<Config>::dryrun(DryRun<typename Configuration::value_type, SignType>& func)
{
    typedef typename Configuration::value_type Vertex;
    for (unsigned int k = 0; k < len; ++k)
    {
        func.template onSector<UP>(k, t_M, k, t_M);
        func.template onSector<DOWN>(k, t_M, k, t_M);
    }
}

/**
A class for measuring the Total Double Occupancy at a certain hardcoded time
*/
template <class Config>
class TotalDoubleOccupancy : public Network_Cache<Config, typename Config::SignType>
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef GFRetVal ObservableType;///< the ParticleNumber is essentially a floating Point type, but during calculations it can be complex
    /**
    The Constructor for the Average Order
    */
    TotalDoubleOccupancy(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, typename Config::SignType>(n, "TotalDoubleOccupancy"), t_M(params.t_exp), len(params.N)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    This determines the Average Order for the given Configuration
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const double& t_M;
    const uint32_t& len;
};

template <class Config>
void TotalDoubleOccupancy<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    typedef typename Configuration::value_type Vertex;
    Vertex v0u(0, 0, UP);
    Vertex v0d(0, 0, DOWN);
    genericTwoParticleGreensfunction_dry<Configuration>(func, v0u, v0u, v0d, v0d);
    for (unsigned int k = 1; k < len; ++k)
    {
      Vertex vu(k, 0, UP);
      Vertex vd(k, 0, DOWN);
      genericTwoParticleGreensfunction_dry<Configuration>(func, vu, vu, vd, vd);
    }
}

template <class Config>
void TotalDoubleOccupancy<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& func)
{
        typename Configuration::value_type v0u(0, 0, UP);
      typename Configuration::value_type v0d(0, 0, DOWN);
    ObservableType obs(genericTwoParticleGreensfunction<Configuration>(func, v0u, v0u, v0d, v0d));
    for (unsigned int k = 1; k < len; ++k)
    {
      typename Configuration::value_type vu(k, 0, UP);
      typename Configuration::value_type vd(k, 0, DOWN);
      obs += genericTwoParticleGreensfunction<Configuration>(func, vu, vu, vd, vd);
    }
//        obs += func.template onSector<UP>(k, t_M, k, t_M) * func.template onSector<DOWN>(k, t_M, k, t_M);
    obs *= configuration.phase;
    this->add_bin(obs);
    return;
}

/**
A class for measuring the time dependency of the kinetic energy as defined in the Hubbard model
*/
template <class Config>
class KineticEnergy : public Network_Cache<Config, std::valarray<typename Config::SignType> >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> ObservableType;///< Kinetic Energy is a function-like observable in realtime evolution
    /**
    The Constructor for the Average Order
    */
    KineticEnergy(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "KineticEnergy"), len(params.N), functionpoints(params.functionpoints), t(params.t), delta_s(params.delta_s)
    {}
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    This determines the Average Order for the given Configuration
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double& t;
    const double delta_s;
};

template <class Config>
void KineticEnergy<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    typedef typename Configuration::value_type Vertex;

    for (unsigned int j = 0; j < functionpoints; ++j)
    {
        const FPType s = j * delta_s;
        for (unsigned int k = 0; k < len; ++k)
        {
            func.template onSector<UP>(k, s, (k + 1)%len, s);
            func.template onSector<UP>(k, s, (len + k - 1)%len, s);

            func.template onSector<DOWN>(k, s, (k + 1)%len, s);
            func.template onSector<DOWN>(k, s, (len + k - 1)%len, s);
        }
    }
}

template <class Config>
void KineticEnergy<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    std::valarray<GFRetVal> vals(functionpoints);
    for (unsigned int j = 0; j < functionpoints; ++j)
    {
        FPType s = j * delta_s;
        GFRetVal obs = 0;
        for (unsigned int k = 0; k < len; ++k)
        {
            obs += dowick.template onSector<UP>(k, s, (k + 1)%len, s);
            obs += dowick.template onSector<UP>(k, s, (len + k - 1)%len, s);

            obs += dowick.template onSector<DOWN>(k, s, (k + 1)%len, s);
            obs += dowick.template onSector<DOWN>(k, s, (len + k - 1)%len, s);
        }
        obs *= static_cast<FPType>(-t) * configuration.phase;
        vals[j] = obs;
    }
    //add to measurement
    this->add_bin(vals);
    return;
}

/**
A class for measuring the time dependency of the Magnetization
*/
template <class Config>
class Magnetization : public Network_Cache<Config, std::valarray<typename Config::SignType> >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> ObservableType;///< the Magnetization is a function of the time
    /**
    The Constructor for the Magnetization
    */
    Magnetization(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "Magnetization"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Magnetization for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void Magnetization<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    typedef typename Configuration::value_type Vertex;

    for (unsigned int j = 0; j < functionpoints; ++j)
    {
        FPType s = j * delta_s;
        for (unsigned int k = 0; k < len; ++k)
        {
            func.template onSector<UP>(k, s, k, s);
            func.template onSector<DOWN>(k, s, k, s);
        }
    }
}

template <class Config>
void Magnetization<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    std::valarray<GFRetVal> vals(functionpoints);
    for (unsigned int j = 0; j < functionpoints; ++j)
    {
        FPType s = j * delta_s;
        GFRetVal obs = 0;
        for (unsigned int k = 0; k < len; ++k)
        {
            obs += dowick.template onSector<UP>(k, s, k, s) - dowick.template onSector<DOWN>(k, s, k, s);
        }
        obs *= configuration.phase;
        vals[j] = obs;
    }
    //add to measurement
    this->add_bin(vals);
    return;
}

/**
A class for measuring the time dependency of the Eta-Pairing
*/
template <class Config>
class EtaPairing : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;//Eta-Pairing in k-space can be complex
    typedef std::valarray<Function> ObservableType;///< The Eta-Pairing is a time-dependent
    /**
    The Constructor for the Eta-Pairing
    */
    EtaPairing(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "EtaPairing"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Eta-Pairing for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void EtaPairing<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (int q = 0; q < static_cast<int>(len); ++q)//for every k-space-value q
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every timeslice
        {
            const FPType s = j * delta_s;//the realtime
            for (int a = 0; a < static_cast<int>(len); ++a)
                for (int d = 0; d < static_cast<int>(len); ++d)
                {
                    func.template onSector<UP>(a, s, d, s);
                    func.template onSector<DOWN>(a, s, d, s);
                }
        }
    }
    return;
}

template <class Config>
void EtaPairing<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (int q = 0; q < static_cast<int>(len); ++q)//for every k-space-value q
    {
        func[q].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)//for every timeslice
        {
            const FPType s = j * delta_s;//the realtime
            std::complex<FPType> t1 = 0;
            FPType lenf = static_cast<FPType>(len);
            for (int a = 0; a < static_cast<int>(len); ++a)
                for (int d = 0; d < static_cast<int>(len); ++d)
                {
                    std::complex<FPType> pref = std::exp(std::complex<FPType>(0.0, static_cast<FPType>(TWOPI*q)/lenf*(d-a)) );
                    t1 += pref * dowick.template onSector<UP>(a, s, d, s) * dowick.template onSector<DOWN>(a, s, d, s);
                }
            //add to measurement
            func[q][j] = t1;
        }
        func[q] *= configuration.phase;
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the Eta-Pairing
*/
template <class Config>
class EtaPairing_Real : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< The Eta-Pairing is a time-dependent observable
    /**
    The Constructor for the Eta-Pairing
    */
    EtaPairing_Real(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "EtaPairing_Real"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Eta-Pairing for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void EtaPairing_Real<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (int q = 0; q < static_cast<int>(len); ++q)//for every lattice site q
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every timeslice
        {
            const FPType s = j * delta_s;//the realtime
            func.template onSector<UP>(q, s, 0, s);
            func.template onSector<DOWN>(q, s, 0, s);
        }
    }
    return;
}

template <class Config>
void EtaPairing_Real<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (uint32_t q = 0; q < len; ++q)
    {
        func[q].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)//for every timeslice
        {
            const FPType s = j * delta_s;//the realtime
            GFRetVal t1 = dowick.template onSector<UP>(q, s, 0, s) * dowick.template onSector<DOWN>(q, s, 0, s);
            //add to measurement
            func[q][j] = t1;
        }
        func[q] *= configuration.phase;
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the Spin-Spin-Correlation
*/
template <class Config>
class SpinSpinCorrelation : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< Spin-Spin-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the Spin-Spin-Correlation
    */
    SpinSpinCorrelation(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SpinSpinCorrelation"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Spin-Spin-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void SpinSpinCorrelation<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, DOWN);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, DOWN);
        }
    }
    return;
}

template <class Config>
void SpinSpinCorrelation<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
            GFRetVal t1 = densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, UP)
                          - densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, DOWN)
                          - densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, UP)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, DOWN);
            //add to measurement
            func[k][j] = t1 * configuration.phase/ static_cast<FPType>(4.0);
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the correlated part of the Spin-Spin-Correlation
*/
template <class Config>
class SpinSpinCorrelatedPart : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< Spin-Spin-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the correlated part Spin-Spin-Correlation
    */
    inline SpinSpinCorrelatedPart(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SpinSpinCorrelatedPart"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Spin-Spin-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void SpinSpinCorrelatedPart<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, DOWN);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, DOWN);
            func.template onSector<UP>(0, s, 0, s);
            func.template onSector<DOWN>(k, s, k, s);
        }
    }
    return;
}

template <class Config>
void SpinSpinCorrelatedPart<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
            GFRetVal au = dowick.template onSector<UP>(0, s, 0, s);
            GFRetVal ad = dowick.template onSector<DOWN>(0, s, 0, s);
            GFRetVal bu = dowick.template onSector<UP>(k, s, k, s);
            GFRetVal bd = dowick.template onSector<DOWN>(k, s, k, s);
            GFRetVal t1 = densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, UP)
                          - densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, DOWN)
                          - densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, UP)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, DOWN)
                          -au*bu + au*bd + ad * bu - ad * bd;
            //add to measurement
            func[k][j] = t1 * configuration.phase/ static_cast<FPType>(4.0);
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the Charge-Charge-Correlation
*/
template <class Config>
class ChargeChargeCorrelation : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< Charge-Charge-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the Charge-Charge-Correlation
    */
    ChargeChargeCorrelation(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "ChargeChargeCorrelation"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Charge-Charge-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void ChargeChargeCorrelation<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, DOWN);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, DOWN);
        }
    }
    return;
}

template <class Config>
void ChargeChargeCorrelation<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
            GFRetVal t1 = densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, UP)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, DOWN)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, UP)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, DOWN);
            //add to measurement
            func[k][j] = t1 * configuration.phase;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the correlated part of the Charge-Charge-Correlation
*/
template <class Config>
class ChargeChargeCorrelatedPart : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< Charge-Charge-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the correlated part of the Charge-Charge-Correlation
    */
    ChargeChargeCorrelatedPart(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "ChargeChargeCorrelatedPart"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the correlated part of the Charge-Charge-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void ChargeChargeCorrelatedPart<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, UP, k, s, DOWN);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, UP);
            densitydensityCorrelation_dry<Configuration>(func, 0, s, DOWN, k, s, DOWN);
            func.template onSector<UP>(0, s, 0, s);
            func.template onSector<DOWN>(k, s, k, s);
        }
    }
    return;
}

template <class Config>
void ChargeChargeCorrelatedPart<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
            GFRetVal au = dowick.template onSector<UP>(0, s, 0, s);
            GFRetVal ad = dowick.template onSector<DOWN>(0, s, 0, s);
            GFRetVal bu = dowick.template onSector<UP>(k, s, k, s);
            GFRetVal bd = dowick.template onSector<DOWN>(k, s, k, s);
            GFRetVal t1 = densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, UP)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, UP, k, s, DOWN)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, UP)
                          + densitydensityCorrelation<Configuration>(dowick, 0, s, DOWN, k, s, DOWN)
                          - au*bu
                          - ad*bu
                          - au*bd
                          - ad*bd;
            //add to measurement
            func[k][j] = t1 * configuration.phase;
        }
    }
    this->add_bin(func);
    return;
}


/**
A class for measuring the time dependency of the k-space resolved particle density. the dependency on the spin is summed out.
It measures \sum_r <c_r^\dagger c_0 >
*/
template <class Config>
class kSpaceDensity : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;
    typedef std::valarray<Function> ObservableType;///< Charge-Charge-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the Charge-Charge-Correlation
    */
    kSpaceDensity(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "kSpaceDensity"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Charge-Charge-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void kSpaceDensity<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int r = 0; r < len; ++r)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
            func.template onSector<UP>(r, s, 0, s);
            func.template onSector<DOWN>(r, s, 0, s);
        }
    }
    return;
}

template <class Config>
void kSpaceDensity<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
            std::complex<FPType> t1 = 0;
            for (unsigned int r = 0; r < len; ++r)
            {
                std::complex<FPType> t2 = dowick.template onSector<UP>(r, s, 0, s) + dowick.template onSector<DOWN>(r, s, 0, s);
                t2 = t2 * std::exp(std::complex<FPType>(0.0, static_cast<FPType>(2.0 * M_PIl/len) * k * r) );
                t1 += t2;
            }
            //add to measurement
            func[k][j] = t1 * configuration.phase;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the LocalDensityVariance
*/
template <class Config>
class LocalDensityVariance : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< LocalDensityVariance is a spatially resolved time-dependent observable
    /**
    The Constructor for the LocalDensityVariance
    */
    LocalDensityVariance(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "LocalDensityVariance"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the LocalDensityVariance for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void LocalDensityVariance<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
            func.template onSector<UP>(k, s, k, s);
            func.template onSector<DOWN>(k, s, k, s);
            densitydensityCorrelation_dry<Configuration>(func, k, s, UP, k, s, DOWN);
        }
    }
    return;
}

template <class Config>
void LocalDensityVariance<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;

            GFRetVal nup = dowick.template onSector<UP>(k, s, k, s);
            GFRetVal ndown = dowick.template onSector<UP>(k, s, k, s);
            GFRetVal nupndown = densitydensityCorrelation<Configuration>(dowick, k, s, UP, k, s, DOWN);
            //add to measurement
            func[k][j] = configuration.phase * (nup*(static_cast<FPType>(1.0) - nup) + ndown*(static_cast<FPType>(1.0) - ndown) + static_cast<FPType>(2.0) * (nupndown - nup*ndown));
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the Greensfunctions, thus <c_0(0)^\dagger c_r(s)>,
WE MEASURE THE UP-SECTOR!
*/
template <class Config>
class Greensfunction : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< The Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    Greensfunction(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType> (n, "Greensfunction"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void Greensfunction<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
/*    	  typename Configuration::value_type v1(k, s, UP);
	  typename Configuration::value_type v2(0, 0, DOWN);
	    func(v1, v2);*/
            func.template onSector<UP>(k, s, 0, 0);
//            func.template onSector<DOWN>(k, s, 0, 0);
	    
        }
    }
    return;
}

template <class Config>
void Greensfunction<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
/*	    
       	  typename Configuration::value_type v1(k, s, UP);
	  typename Configuration::value_type v2(0, 0, DOWN);
	  GFRetVal t1 = dowick(v1, v2);*/
            GFRetVal t1 = dowick.template onSector<UP>(k, s, 0, 0)
            //+dowick.template onSector<DOWN>(k, s, 0, 0)
                          ;
            //add to measurement
            func[k][j] = configuration.phase * t1;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the Diagonal Green's function <c_ks(tau)^\dagger c_ks(0)>,
We employ Time-reversal symmetry to enhance the measurement. In theory it is always a real quantity
*/
template <class Config>
class DiagonalGreensfunction_kspace : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;
    typedef std::valarray<Function> ObservableType;///< The Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    DiagonalGreensfunction_kspace(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType> (n, "Greensfunction"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void DiagonalGreensfunction_kspace<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int r = 0; r < len; ++r)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
	 func.template onSector<UP>(r, s, 0, 0);
	 func.template onSector<DOWN>(r, s, 0, 0);
        }
    }
    return;
}

template <class Config>
void DiagonalGreensfunction_kspace<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    FPType invlen = static_cast<FPType>(1.0)/len;
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
	    std::complex<FPType> sum = 0;
            for(uint r = 0; r < len; ++r)
            {
                std::complex<FPType> a = std::exp(std::complex<FPType>(0.0, -static_cast<FPType>(TWOPI*k*r)*invlen));
		GFRetVal t1 = dowick.template onSector<UP>(r, s, 0, 0) + dowick.template onSector<DOWN>((len-r)%len, s, 0, 0); 
                sum += a * t1;
            }
            func[k][j] = 0.5*configuration.phase * sum /* invlen*/;//not sure about that invlen here...
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the TRI Diagonal Green's function <c_ks(tau)^\dagger c_ks(0)>,
s = +-
*/
template <class Config, int sign>
class TRI_Greensfunction_kspace : public Network_Cache<Config, std::valarray<std::valarray<
//typename Config::Configuration::FPType
typename Config::SignType
> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;//G++ and G-- are real
    typedef std::valarray<Function> ObservableType;///< The Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    TRI_Greensfunction_kspace(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType> (n, "TRI_Greensfunction_kspace"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config, int sign>
void TRI_Greensfunction_kspace<Config, sign>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    typename Configuration::value_type v2(0, 0, DOWN);
    typename Configuration::value_type v4(0, 0, UP);
    for (unsigned int j = 0; j < functionpoints; ++j)
    {
            const FPType s = j * delta_s;
            for(uint r = 0; r < len; ++r)
            {
		typename Configuration::value_type v1(r, s, UP);
		typename Configuration::value_type v3(r, s, DOWN);
		func.template onSector<UP>(r, s, 0, 0);
//		func.template onSector<DOWN>(r, s, 0, 0);
		func(v1, v2);
//		func(v3, v4);
            }
    }
    return;
}

template <class Config, int sign>
void TRI_Greensfunction_kspace<Config, sign>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    FPType invlen = static_cast<FPType>(1.0)/len;
    typename Configuration::value_type v2(0, 0, DOWN);
    typename Configuration::value_type v4(0, 0, UP);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
	    GFRetVal sum = 0;
            for(uint r = 0; r < len; ++r)
            {
                std::complex<FPType> a = std::exp(std::complex<FPType>(0.0, -static_cast<FPType>(TWOPI*k*r)*invlen));
		typename Configuration::value_type v1(r, s, UP);
		typename Configuration::value_type v3(r, s, DOWN);
//		std::complex<FPType> t1 = dowick.template onSector<UP>(r, s, 0, 0) + dowick.template onSector<DOWN>(r, s, 0, 0) 
//		             + std::complex<FPType>(0.0, sign)*( dowick(v1, v2) - dowick(v3, v4));
//the next line is equivalent down to the Monte-Carlo level it seems...
//		std::complex<FPType> t1 = dowick.template onSector<UP>(r, s, 0, 0) + std::complex<FPType>(0.0, sign)*dowick(v1,v2);
		GFRetVal gd = dowick.template onSector<UP>(r, s, 0, 0);//measure diagonal
		GFRetVal go = dowick(v1, v2);//measure offdiagonal
		sum += real(a)*gd - sign * imag(a) * go;
            }
            func[k][j] = configuration.phase * sum;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the Diagonal Green's function <c_ks(tau)^\dagger c_k(-s)(0)>,
We employ Time-reversal symmetry to enhance the measurement. 
*/
template <class Config>
class OffDiagonalGreensfunction_kspace : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;
    typedef std::valarray<Function> ObservableType;///< The Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    OffDiagonalGreensfunction_kspace(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType> (n, "Greensfunction"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void OffDiagonalGreensfunction_kspace<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
  typename Configuration::value_type v2(0, 0, DOWN);
  typename Configuration::value_type v4(0, 0, UP);
    for (unsigned int r = 0; r < len; ++r)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
    	  typename Configuration::value_type v1(r, s, UP);
	  func(v1, v2);
	  typename Configuration::value_type v3(r, s, DOWN);
	  func(v3, v4);
        }
    }
    return;
}

template <class Config>
void OffDiagonalGreensfunction_kspace<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    FPType invlen = static_cast<FPType>(1.0)/len;
    typename Configuration::value_type v2(0, 0, DOWN);
    typename Configuration::value_type v4(0, 0, UP);
    const FPType fac = TWOPI * invlen;
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
	    std::complex<FPType> sum = 0;
            for(uint r = 0; r < len; ++r)
            {
                std::complex<FPType> pref = std::exp(std::complex<FPType>(0.0, -fac * (k*r)));
       	        typename Configuration::value_type v1(r, s, UP);
		typename Configuration::value_type v3((len - r)%len, s, DOWN);
	        GFRetVal t1 = dowick(v1, v2) + dowick(v3, v4);
                sum += pref * t1;
            }
            func[k][j] = 0.5*configuration.phase * sum;//normalization not necessary since the 1/N factor is already in the definition of G_0
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the Greensfunctions, thus <c_0(0)^\dagger c_r(s)>.
This Greensfunction only has sense for imaginary time measurements, because we use David's smoothing trick.
WE MEASURE THE DOWN-SECTOR!
*/
template <class Config>
class SmoothImaginaryGreensfunction : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< The Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    SmoothImaginaryGreensfunction(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SmoothImaginaryGreensfunction"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s), beta(params.beta)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
    const double beta;
};

template <class Config>
void SmoothImaginaryGreensfunction<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            for (unsigned int i = 0; i < trunc(beta/delta_s); ++i)
            {
//            func.template onSector<UP>(k, s, 0, 0);
                func.template onSector<DOWN>(0, static_cast<FPType>(i+j)*delta_s, k, i * delta_s);
            }
        }
    }
    return;
}

template <class Config>
void SmoothImaginaryGreensfunction<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            GFRetVal sum = 0;
            for (unsigned int i = 0; i < trunc(beta/delta_s); ++i)
            {
                sum += /*dowick.template onSector<UP>(k, s, 0, 0)
            +*/ dowick.template onSector<DOWN>(0, static_cast<FPType>(i+j)*delta_s, k, i * delta_s) * delta_s
                    ;
            }
            //add to measurement
            func[k][j] = configuration.phase * sum/beta;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the Greensfunctions, thus <c_0(0)^\dagger c_r(s)>.
This Greensfunction only has sense for imaginary time measurements, because we use David's smoothing trick.
Here we average over both Spin sectors
*/
template <class Config>
class SmoothImaginaryGreensfunction_averaged : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< The Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    SmoothImaginaryGreensfunction_averaged(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SmoothImaginaryGreensfunction_averaged"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s), beta(params.beta)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
    const double beta;
};

template <class Config>
void SmoothImaginaryGreensfunction_averaged<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            for (unsigned int i = 0; i < trunc(beta/delta_s); ++i)
            {
                func.template onSector<DOWN>(0, static_cast<FPType>(i+j)*delta_s, k, i * delta_s);
		func.template onSector<UP>(0, static_cast<FPType>(i+j)*delta_s, k, i * delta_s);
            }
        }
    }
    return;
}

template <class Config>
void SmoothImaginaryGreensfunction_averaged<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            GFRetVal sum = 0;
            for (unsigned int i = 0; i < trunc(beta/delta_s); ++i)
            {
                sum += (dowick.template onSector<DOWN>(0, static_cast<FPType>(i+j)*delta_s, k, i * delta_s) + dowick.template onSector<UP>(0, static_cast<FPType>(i+j)*delta_s, k, i * delta_s)) * delta_s;
            }
            //add to measurement
            func[k][j] = configuration.phase * sum/beta/2.0;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the Matsubarafrequency Greensfunction thus <c_0(0)^\dagger c_r(i omega)>,
WE MEASURE THE UP-SECTOR!
*/
template <class Config, class Greensfunction>
class MatsubaraFrequencyGreensfunction : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;
    typedef std::valarray<Function> ObservableType;///< The Matsubarafrequency dependent Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    MatsubaraFrequencyGreensfunction(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "MatsubaraFrequencyGreensfunction"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s), beta(params.beta)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&) MTL_CONST_FUNCTION;
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const FPType delta_s;
    const FPType beta;
};

template <class Config, class Greensfunction>
void MatsubaraFrequencyGreensfunction<Config, Greensfunction>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>&)
{
}

template <class Config, class Greensfunction>
void MatsubaraFrequencyGreensfunction<Config, Greensfunction>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>&)
{
    ObservableType func(len);
    typename Greensfunction::Vertex v;
    v.spin = UP;
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int r = 0; r < configuration.size(); ++r)
            for (unsigned int s = 0; s < configuration.size(); ++s)
            {
                const FPType deltars = configuration[r].tau - configuration[s].tau;
                std::complex<FPType> inc = std::exp(std::complex<FPType>(0.0, deltars * M_PI / beta));
                const std::complex<FPType> fac = inc * inc;
                inc *= configuration.matcont(r, s, UP, UP);
                for (unsigned int n = 0; n < functionpoints; ++n)
                {
                    func[k][n] += inc;
                    inc *= fac;
                }
            }
        for (unsigned int n = 0; n < functionpoints; ++n)
        {
            FPType omegan = M_PI/beta*static_cast<FPType>(2*n+1);
            std::complex<FPType> gomegan = Greensfunction::gomega(omegan, v, v);
            func[k][n] = configuration.phase * gomegan * (static_cast<FPType>(1.0) - gomegan*func[k][n]/beta);
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the averaged Matsubarafrequency Greensfunction thus sum_\sigma <c_0(0)^\dagger c_r(i omega)>,
*/
template <class Config, class Greensfunction>
class MatsubaraFrequencyGreensfunction_averaged : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;
    typedef std::valarray<Function> ObservableType;///< The Matsubarafrequency dependent Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    MatsubaraFrequencyGreensfunction_averaged(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "MatsubaraFrequencyGreensfunction_averaged"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s), beta(params.beta)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&) MTL_CONST_FUNCTION;
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const FPType delta_s;
    const FPType beta;
};

template <class Config, class Greensfunction>
void MatsubaraFrequencyGreensfunction_averaged<Config, Greensfunction>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>&)
{
}

template <class Config, class Greensfunction>
void MatsubaraFrequencyGreensfunction_averaged<Config, Greensfunction>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>&)
{
    ObservableType func(len);
    typename Greensfunction::Vertex vup;
    typename Greensfunction::Vertex vdown;
    vup.spin = UP;
    vdown.spin = DOWN;
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
	std::complex<FPType> gup[functionpoints];
	std::complex<FPType> gdown[functionpoints];
	for (unsigned int n = 0; n < functionpoints; ++n)
        {
            FPType omegan = M_PI/beta*static_cast<FPType>(2*n+1);
            std::complex<FPType> gomeganup = Greensfunction::gomega(omegan, vup, vup);
	    std::complex<FPType> gomegandown = Greensfunction::gomega(omegan, vdown, vdown);
	    gup[n] = gomeganup;
	    gdown[n] = gomegandown;
        }
	
	
        for (unsigned int r = 0; r < configuration.size(); ++r)
            for (unsigned int s = 0; s < configuration.size(); ++s)
            {
                const FPType deltars = configuration[r].tau - configuration[s].tau;
                std::complex<FPType> inc = std::exp(std::complex<FPType>(0.0, deltars * M_PI / beta));
                const std::complex<FPType> fac = inc * inc;
//                inc *= configuration.matcont.mat(2*r, 2*s);
                for (unsigned int n = 0; n < functionpoints; ++n)
                {
                    func[k][n] += (gup[n]*gup[n]*configuration.matcont.mat(2*r, 2*s) + gdown[n]*gdown[n]*configuration.matcont.mat(2*r+1, 2*s+1))*inc;
                    inc *= fac;
                }
            }
        for (unsigned int n = 0; n < functionpoints; ++n)
        {
            func[k][n] = configuration.phase * (gup[n] + gdown[n] - func[k][n]/beta);
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the derivative with respect to omega of the Matsubarafrequency Greensfunction thus d/domega <c_0(0)^\dagger c_r(i omega)>,
WE MEASURE THE UP-SECTOR!
*/
template <class Config, class Greensfunction>
class MatsubaraFrequencyGreensfunctionDerivative : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;
    typedef std::valarray<Function> ObservableType;///< The Matsubarafrequency dependent Greensfunction contains an array of functions
    /**
    The Constructor for the Greensfunction
    */
    MatsubaraFrequencyGreensfunctionDerivative(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "MatsubaraFrequencyGreensfunctionDerivative"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s), beta(params.beta)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&) MTL_CONST_FUNCTION;
    /**
    Here we evaluate for a given order the values of all Greensfunctions
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const FPType delta_s;
    const FPType beta;
};

template <class Config, class Greensfunction>
void MatsubaraFrequencyGreensfunctionDerivative<Config, Greensfunction>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>&)
{
}

template <class Config, class Greensfunction>
void MatsubaraFrequencyGreensfunctionDerivative<Config, Greensfunction>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>&)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        Function temp(functionpoints);
        /*                for (unsigned int n = 0; n < functionpoints; ++n)
                        {
                            const FPType omegan = M_PI/beta*static_cast<FPType>(2*n+1);
                            const std::complex<FPType> gomegan = CRAPPY_Make_compile_Helper<Greensfunction, Configuration>::measure(omegan, beta, v, w, ed, mu);
                            for(unsigned int r = 0; r < configuration.size(); ++r)
                            {
                              for(unsigned int s = 0; s < configuration.size(); ++s)
                                  func[k][n] += std::exp(std::complex<FPType>(0.0, omegan * (configuration[r].tau - configuration[s].tau))) * configuration.up.inverse(r,s);
                            }
                            func[k][n] /= beta;
                            func[k][n] *= gomegan;
                            func[k][n] = configuration.phase * gomegan * (1.0 - func[k][n]);
                        }*/

        for (unsigned int r = 0; r < configuration.size(); ++r)
            for (unsigned int s = 0; s < configuration.size(); ++s)
            {
                const FPType deltars = configuration[r].tau - configuration[s].tau;
                std::complex<FPType> inc1 = std::exp(std::complex<FPType>(0.0, deltars * M_PI / beta));
                const std::complex<FPType> fac = inc1 * inc1;
                inc1 *= configuration.matcont.up(r, s);
                std::complex<FPType> inc2 = inc1 * std::complex<FPType>(0.0, deltars);
                for (unsigned int n = 0; n < functionpoints; ++n)
                {
                    func[k][n] += inc1;//func[k][n] is the sum, that is the same as in the plain Matsubara Greensfunction
                    temp[n] += inc2;
                    inc1 *= fac;
                    inc2 *= fac;
                }
            }
        for (unsigned int n = 0; n < functionpoints; ++n)
        {
            FPType omegan = M_PI/beta*static_cast<FPType>(2*n+1);
            const FPType tiny = 0.00000001;//same as in the greensfunctions
            const FPType h = pow(tiny, 1.0/3.0) * 0.707 * omegan;//an optimal choice of h for the symmetric derivative derived for a function that behaves as 1/x... see NR 5.7
            std::complex<FPType> gomegan = Greensfunction::gomega(omegan, UP);
            std::complex<FPType> gomeganplush = Greensfunction::gomega(omegan + h, UP);
            std::complex<FPType> gomeganminush = Greensfunction::gomega(omegan - h, UP);
            std::complex<FPType> derivative = (gomeganplush - gomeganminush)/(static_cast<FPType>(2.0)*h);
            func[k][n] = configuration.phase * (derivative - static_cast<FPType>(2.0)/beta * gomegan * derivative * func[k][n] - static_cast<FPType>(1.0)/beta* gomegan * gomegan * temp[n]);
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the Density-Density-StructureFactor
*/
template <class Config>
class DensityDensityStructureFactor : public Network_Cache<Config, std::valarray<std::valarray<std::complex<typename Config::Configuration::FPType> > > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<std::complex<FPType> > Function;//the density-density structure factor is determined by a Fourier-transform. Therefore it can be a complex type
    typedef std::valarray<Function> ObservableType;///< the density-density-structure-factor is a k-space resolved time-dependent observables
    /**
    The Constructor for the DensityDensityStructureFactor
    */
    DensityDensityStructureFactor(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "DensityDensityStructureFactor"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    This determines the DensityDensityStructureFactor for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void DensityDensityStructureFactor<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)
    {
        for (unsigned int a = 0; a < len; ++a)
        {
            for (unsigned int j = 0; j < functionpoints; ++j)
            {
                const FPType s = j * delta_s;
                densitydensityCorrelation_dry<Configuration>(func, a, s, UP, k, s, UP);
                densitydensityCorrelation_dry<Configuration>(func, a, s, UP, k, s, DOWN);
                densitydensityCorrelation_dry<Configuration>(func, a, s, DOWN, k, s, UP);
                densitydensityCorrelation_dry<Configuration>(func, a, s, DOWN, k, s, DOWN);
            }
        }
    }
    return;
}

template <class Config>
void DensityDensityStructureFactor<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
            std::complex<FPType> t1 = 0;
            FPType lenf = static_cast<FPType>(len);
            for (int a = 0; a < static_cast<int>(len); ++a)
                for (int d = 0; d < static_cast<int>(len); ++d)
                {
                    std::complex<FPType> pref = std::exp(std::complex<FPType>(0.0, static_cast<FPType>(TWOPI*k)/lenf*(d-a)) );
                    t1 +=     pref * (densitydensityCorrelation<Configuration>(dowick, a, s, UP, d, s, UP)
                                      + densitydensityCorrelation<Configuration>(dowick, a, s, UP, d, s, DOWN)
                                      + densitydensityCorrelation<Configuration>(dowick, a, s, DOWN, d, s, UP)
                                      + densitydensityCorrelation<Configuration>(dowick, a, s, DOWN, d, s, DOWN));
                }
            //add to measurement
            func[k][j] = t1 * configuration.phase/lenf;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the DoubleOccupancy
*/
template <class Config>
class DoubleOccupancy : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< DoubleOccupancy is a spatially resolved time-dependent observable
    /**
    The Constructor for the DoubleOccupancy
    */
    DoubleOccupancy(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config,ObservableType>(n, "DoubleOccupancy"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the DoubleOccupancy for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void DoubleOccupancy<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    for (unsigned int k = 0; k < len; ++k)//for each site
    {
        for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
        {
            const FPType s = j * delta_s;
            densitydensityCorrelation_dry<Configuration>(func, k, s, DOWN, k, s, UP);
        }
    }
    return;
}

template <class Configuration>
void DoubleOccupancy<Configuration>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
            GFRetVal t1 = densitydensityCorrelation<Configuration>(dowick, k, s, DOWN, k, s, UP);
            //add to measurement
            func[k][j] = t1 * configuration.phase;
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the Conductance in an imaginary time simulation
*/
#include "conductance_grid.h"
template <class Config, class Greensfunction>
class Conductance : public Network_Cache<Config, typename Config::SignType>
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef GFRetVal ObservableType;
    Conductance(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "Conductance"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s), beta(params.beta), alpha_max(sizeof(grid)/sizeof(Point2)), gamma(params.V*params.V*M_PI/params.W)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&) MTL_CONST_FUNCTION;
    /**
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
    const FPType beta;
    const unsigned int alpha_max;
    const FPType gamma;
};

template <class Configuration, class Greensfunction>
void Conductance<Configuration, Greensfunction>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    return;
}

template <class Configuration, class Greensfunction>
void Conductance<Configuration, Greensfunction>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType ret = 0.0;
    for (unsigned int alpha = 0; alpha < alpha_max; ++alpha)
    {
        const FPType omega_alpha = 1.0/(grid[alpha].gp*beta);
        const std::complex<FPType> gomega = Greensfunction::gomega(omega_alpha, UP);
        const FPType tiny = 0.00000001;//same as in the greensfunctions
        const FPType h = pow(tiny, 1.0/3.0) * 0.707 * omega_alpha;//an optimal choice of h for the symmetric derivative, derived for a function that behaves as 1/x... see NR 5.7. the choice of x_c= 0.707*x is reasonable as we always evaluate the function for x > 0
        std::complex<FPType> gomegaplush = Greensfunction::gomega(omega_alpha + h, UP);
        std::complex<FPType> gomegaminush = Greensfunction::gomega(omega_alpha - h, UP);
        std::complex<FPType> derivative = (gomegaplush - gomegaminush)/(static_cast<FPType>(2.0)*h);
        std::complex<FPType> sum1 = 0.0;
        std::complex<FPType> sum2 = 0.0;
        for (unsigned int r = 0; r < configuration.size(); ++r)
        {
            for (unsigned int s = 0; s < configuration.size(); ++s)
            {
                std::complex<FPType> expfac = std::exp(std::complex<FPType>(0.0, omega_alpha * (configuration[r].tau - configuration[s].tau)));
                GFRetVal mat = configuration.matcont.up(r,s);
                sum1 += mat*expfac;
                sum2 += mat*expfac * std::complex<FPType>(0.0, (configuration[r].tau - configuration[s].tau));
            }
        }
        std::complex<FPType> factor = gomega / beta;
        std::complex<FPType> dga = derivative - static_cast<FPType>(2.0) * factor * derivative * sum1 - factor *gomega * sum2;
        ret += imag(dga) * grid[alpha].weight;
    }
    this->add_bin(-configuration.phase * 2.0*gamma*ret/beta);
}

/**
A class for measuring the time dependency of the correlated part of the Spin-Spin-Correlation
*/
template <class Config>
class SpinSpinCorrelatedPart_Y : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< Spin-Spin-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the correlated part Spin-Spin-Correlation
    */
    inline SpinSpinCorrelatedPart_Y(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SpinSpinCorrelatedPart_Y"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Spin-Spin-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void SpinSpinCorrelatedPart_Y<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    typedef typename Configuration::value_type Vertex;
    for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
    {
        for (unsigned int k = 0; k < len; ++k)//for each site
        {
            const FPType s = j * delta_s;
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, UP), Vertex(0, 0, DOWN));
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, UP), Vertex(0, 0, DOWN));
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, DOWN), Vertex(0, 0, UP));
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, DOWN), Vertex(0, 0, UP));
            func(Vertex(k, s, UP), Vertex(k, s, DOWN));
            func(Vertex(k, s, DOWN), Vertex(k, s, UP));
            func(Vertex(0, 0, UP), Vertex(0, 0, DOWN));
            func(Vertex(0, 0, DOWN), Vertex(0, 0, UP));
        }
    }
    return;
}

template <class Config>
void SpinSpinCorrelatedPart_Y<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    typedef typename Configuration::value_type Vertex;
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
	    GFRetVal retval = -(genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, UP), Vertex(0, 0, DOWN))
	    -genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, UP), Vertex(0, 0, DOWN))
	    -genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, DOWN), Vertex(0, 0, UP))
	    +genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, DOWN), Vertex(0, 0, UP)))
	    +/*subtract correlated part, additional minus sign due to i^2*/
	    (dowick(Vertex(k, s, UP), Vertex(k, s, DOWN))-
	    dowick(Vertex(k, s, DOWN), Vertex(k, s, UP)))*
	    (dowick(Vertex(0, 0, UP), Vertex(0, 0, DOWN))-
	    dowick(Vertex(0, 0, DOWN), Vertex(0, 0, UP)));
            //add to measurement
            func[k][j] = retval * configuration.phase/ static_cast<FPType>(4.0);
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the correlated part of the Spin-Spin-Correlation
*/
template <class Config>
class SpinSpinCorrelatedPart_X : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< Spin-Spin-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the correlated part Spin-Spin-Correlation
    */
    inline SpinSpinCorrelatedPart_X(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SpinSpinCorrelatedPart_X"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Spin-Spin-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void SpinSpinCorrelatedPart_X<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    typedef typename Configuration::value_type Vertex;
    for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
    {
        for (unsigned int k = 0; k < len; ++k)//for each site
        {
            const FPType s = j * delta_s;
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, UP), Vertex(0, 0, DOWN));
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, UP), Vertex(0, 0, DOWN));
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, DOWN), Vertex(0, 0, UP));
            genericTwoParticleGreensfunction_dry<Configuration>(func, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, DOWN), Vertex(0, 0, UP));
            func(Vertex(k, s, UP), Vertex(k, s, DOWN));
            func(Vertex(k, s, DOWN), Vertex(k, s, UP));
            func(Vertex(0, 0, UP), Vertex(0, 0, DOWN));
            func(Vertex(0, 0, DOWN), Vertex(0, 0, UP));
        }
    }
    return;
}

template <class Config>
void SpinSpinCorrelatedPart_X<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    typedef typename Configuration::value_type Vertex;
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
/*	    std::cout<<genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, UP), Vertex(0, 0, DOWN))
	    <<" "<<genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, UP), Vertex(0, 0, DOWN))<<" "
	    <<genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, DOWN), Vertex(0, 0, UP))<<" "
	    <<genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, DOWN), Vertex(0, 0, UP))<<std::endl;*/
//simplified measurement in case of spin-diagonal measurement.
	    GFRetVal retval = dowick(Vertex(k,s,UP), Vertex(k,0,UP)) * dowick(Vertex(k,0,DOWN), Vertex(k, s, DOWN)) + 
	    dowick(Vertex(k, s, DOWN), Vertex(k, 0, DOWN)) * dowick(Vertex(k, 0, UP), Vertex(k, s, UP))
	    /*(
	    genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, UP), Vertex(0, 0, DOWN))
	    +genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, UP), Vertex(0, 0, DOWN))
	    +genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, UP), Vertex(k, s, DOWN), Vertex(0, 0, DOWN), Vertex(0, 0, UP))
	    +genericTwoParticleGreensfunction<Configuration>(dowick, Vertex(k, s, DOWN), Vertex(k, s, UP), Vertex(0, 0, DOWN), Vertex(0, 0, UP))
	    )*/
	    //subtract correlated part
/*	    -(dowick(Vertex(k, s, UP), Vertex(k, s, DOWN))+
	    dowick(Vertex(k, s, DOWN), Vertex(k, s, UP)))*
	    (dowick(Vertex(0, 0, UP), Vertex(0, 0, DOWN))+
	    dowick(Vertex(0, 0, DOWN), Vertex(0, 0, UP)))*/;
            //add to measurement
//	    std::cout<<retval<<std::endl;
            func[k][j] = retval * configuration.phase/ static_cast<FPType>(4.0);
        }
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the time dependency of the correlated part of the Spin-Spin-Correlation
*/
template <class Config>
class ImaginarySpinSpinCorrelation_Z : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///< Spin-Spin-Correlations are spatially resolved time-dependent observables
    /**
    The Constructor for the correlated part Spin-Spin-Correlation
    */
    inline ImaginarySpinSpinCorrelation_Z(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "ImaginarySpinSpinCorrelation_Z"), len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Spin-Spin-Correlation for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void ImaginarySpinSpinCorrelation_Z<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    typedef typename Configuration::value_type Vertex;
    for (unsigned int j = 0; j < functionpoints; ++j)//for every time-slice
    {
        for (unsigned int k = 0; k < len; ++k)//for each site
        {
            const FPType s = j * delta_s;
	    densitydensityCorrelation_dry<Configuration>(func, k, s, UP, 0, 0, UP);
	    densitydensityCorrelation_dry<Configuration>(func, k, s, DOWN, 0, 0, DOWN);
	    densitydensityCorrelation_dry<Configuration>(func, k, s, UP, 0, 0, DOWN);
	    densitydensityCorrelation_dry<Configuration>(func, k, s, DOWN, 0, 0, UP);
        }
    }
    return;
}

template <class Config>
void ImaginarySpinSpinCorrelation_Z<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    typedef typename Configuration::value_type Vertex;
    ObservableType func(len);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            const FPType s = j * delta_s;
	    GFRetVal retval = densitydensityCorrelation<Configuration>(dowick, k, s, UP, 0, 0, UP)
	    + densitydensityCorrelation<Configuration>(dowick, k, s, DOWN, 0, 0, DOWN)
	    - densitydensityCorrelation<Configuration>(dowick, k, s, UP, 0, 0, DOWN)
	    - densitydensityCorrelation<Configuration>(dowick, k, s, DOWN, 0, 0, UP);
            //add to measurement
            func[k][j] = retval * configuration.phase/ static_cast<FPType>(4.0);
        }
    }
    this->add_bin(func);
    return;
}

/**
The Spin-Susceptiblity in X direction
*/
template <class Config, class Greensfunction>
class SpinSusceptibility_X : public Network_Cache<Config, typename Greensfunction::GOmegaRetType>
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef typename Greensfunction::GOmegaRetType GOmegaRetType;
    typedef GOmegaRetType ObservableType;
    /**
    */
    inline SpinSusceptibility_X(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SpinSusceptibility_X"), len(params.N), beta(params.beta)
    {
        myf = new GOmegaRetType*[8];
        delta_tau = beta/points;
        const uint tablesize = 2 * points;
        for (uint k = 0; k < 8; ++k)
        {
            myf[k] = new GOmegaRetType[tablesize];
            memset(myf[k], 0, tablesize * sizeof(GOmegaRetType));
        }
        for (uint n = 0; n < points; ++n)
        {
	    FPType omegan = M_PI/beta*(2*n+1);
            for (int s = 0; s < 8; ++s)
            {
                SPINS sigma = (s>>2&1? DOWN: UP);
                SPINS sigma_r = (s>>1&1? DOWN: UP);
                SPINS sigma_s = (s&1? DOWN: UP);
                typename Greensfunction::Vertex v1;
                typename Greensfunction::Vertex v2;
                v1.spin = sigma_r;
                v2.spin = !sigma;
                GOmegaRetType goma = Greensfunction::gomega(omegan, v1, v2);
                GOmegaRetType gomam = Greensfunction::gomega(-omegan, v1, v2);
                v1.spin = sigma;
                v2.spin = sigma_s;
                GOmegaRetType gomb = Greensfunction::gomega(omegan, v1, v2);
                GOmegaRetType gombm = Greensfunction::gomega(-omegan, v1, v2);
                for (uint k = 0; k < points; ++k)
                {
                    FPType tau = k*delta_tau;
                    std::complex<FPType> expt = exp(std::complex<FPType>(0.0, omegan * tau));
                    std::complex<FPType> expmt = std::conj(expt);
                    accessphi(sigma, sigma_r, sigma_s, k) += goma * gomb * expmt + gomam*gombm*expt;
                    accessphi(sigma, sigma_r, sigma_s, points + k) += goma * gomb * expt + gomam*gombm*expmt;
                }
            }
        }
        for (uint s = 0; s < 8; ++s)
            for (uint k = 0; k < tablesize; ++k)
                myf[s][k] /= beta;
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the Spin-Susceptiblity for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (GCC_VERSION > GCC_VER(4,5,0))
    static constexpr uint points = 1000;
#else
    static const uint points = 1000;
#endif
    const uint32_t& len;
    FPType beta;
    FPType delta_tau;
    GOmegaRetType** myf;
    GOmegaRetType& accessphi(SPINS sigma, SPINS sigma_r, SPINS sigma_s, uint tau_idx)
    {
        return myf[
                   (sigma == UP ? 0: 4) +
                   (sigma_r == UP ? 0: 2) +
                   (sigma_s == UP ? 0 : 1)
               ][tau_idx];
    }
    GOmegaRetType contrib(SPINS sigma, SPINS sigmaprime, const Configuration& config, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
    {
        typename Greensfunction::Vertex v1;
        v1.spin = sigmaprime;
        v1.tau = 0;
        typename Greensfunction::Vertex v2;
        v2.spin = !sigmaprime;
        v2.tau = 0;
        GOmegaRetType sum1 = dowick(v1, v2);
        v1.spin = sigma;
        v2.spin = !sigma;
        GFRetVal gf1 = beta * Greensfunction::eval(v1, v2);
        GOmegaRetType suma = 0.0;
        for (uint r = 0; r < config.size(); ++r)
            for (uint s = 0; s < config.size(); ++s)
            {
                FPType taurs = config[r].tau - config[s].tau;
                uint tau_idx = 0;
                if (taurs < 0)
                {
                    tau_idx += points;
                    taurs = -taurs;
                }
                tau_idx += static_cast<uint>(trunc(taurs/delta_tau));
                suma += config.matcont(r, s, UP, UP) * accessphi(sigma, UP, UP, tau_idx);
                suma += config.matcont(r, s, DOWN, UP) * accessphi(sigma, DOWN, UP, tau_idx);
                suma += config.matcont(r, s, UP, DOWN) * accessphi(sigma, UP, DOWN, tau_idx);
                suma += config.matcont(r, s, DOWN, DOWN) * accessphi(sigma, DOWN, DOWN, tau_idx);
            }
        sum1 *= (GOmegaRetType(gf1) - suma);
        GOmegaRetType sum2 = accessphi(sigma, sigmaprime, !sigmaprime, 0);
	uint twosize = static_cast<unsigned int>(2 * config.size());
        typename Configuration::MatConf::MatType psi(1, twosize);
        for (uint k = 0; k < config.size(); ++k)
        {
            v1.tau = config[k].tau;
            v2.tau = 0;
            v1.spin = UP;
            v2.spin = !sigmaprime;
            psi(0,2*k) = Greensfunction::eval(v1, v2);
            v1.spin = DOWN;
            psi(0, 2*k + 1) = Greensfunction::eval(v1, v2);
        }
        config.matcont.multiplyVectorbyConfiguration_right(psi);
//        psi = psi * config.matcont.mat;
        GOmegaRetType sum3 = 0.0;
        for (uint k = 0; k < config.size(); ++k)
        {
            FPType tau = config[k].tau;
            uint tau_idx = points + static_cast<uint>(trunc(tau/delta_tau));/*see formula. that way we take the negative sign into consideration*/
            sum3 += psi(0, 2*k) * accessphi(sigma, UP, sigmaprime, tau_idx);
            sum3 += psi(0, 2*k + 1) * accessphi(sigma, DOWN, sigmaprime, tau_idx);
        }
        GOmegaRetType sum4 = 0.0;
        typename Configuration::MatConf::MatType chi(twosize, 1);
        for (uint k = 0; k < static_cast<uint>(config.size()); ++k)
        {
            v1.tau = 0;
            v2.tau = config[k].tau;
            v2.spin = UP;
            v1.spin = sigmaprime;
            chi(2*k, 0) = Greensfunction::eval(v1, v2);
            v2.spin = DOWN;
            chi(2*k + 1, 0) = Greensfunction::eval(v1, v2);
        }
        config.matcont.multiplyVectorbyConfiguration_left(chi);
//        chi = config.matcont.mat * chi;
        for (uint k = 0; k < static_cast<uint>(config.size()); ++k)
        {
            FPType tau = config[k].tau;
            uint tau_idx = static_cast<uint>(trunc(tau/delta_tau));//tau should be positive being straight from a vertex
            sum4 += chi(2*k, 0) * accessphi(sigma, !sigmaprime, UP, tau_idx);
            sum4 += chi(2*k + 1, 0) * accessphi(sigma, !sigmaprime, DOWN, tau_idx);
        }
        GOmegaRetType sum5 = 0;//we can reuse psi and chi
        for (uint r = 0; r < static_cast<uint>(config.size()); ++r)
            for (uint s = 0; s < static_cast<uint>(config.size()); ++s)
            {
                FPType taurs = config[s].tau - config[r].tau;
                uint tau_idx = 0;
                if (taurs < 0)
                {
                    tau_idx += points;
                    taurs = -taurs;
                }
                tau_idx += static_cast<uint>(trunc(taurs/delta_tau));
                sum5 += psi(0, 2*r) * chi(2*s, 0)*accessphi(sigma, UP, UP, tau_idx);
                sum5 += psi(0, 2*r+1) * chi(2*s, 0)*accessphi(sigma, DOWN, UP, tau_idx);
                sum5 += psi(0, 2*r) * chi(2*s+1, 0)*accessphi(sigma, UP, DOWN, tau_idx);
                sum5 += psi(0, 2*r+1) * chi(2*s+1, 0)*accessphi(sigma, DOWN, DOWN, tau_idx);
            }
        return sum1 - sum2 /*+*/- sum3 /*+*/- sum4 - sum5;
    }
    void contrib_dry(SPINS sigma, SPINS sigmaprime, DryRun<typename Configuration::value_type, GFRetVal>& func)
    {
        typename Greensfunction::Vertex v1;
        v1.spin = sigmaprime;
        v1.tau = 0;
        typename Greensfunction::Vertex v2;
        v2.spin = !sigmaprime;
        v2.tau = 0;
        func(v1, v2);
    }
};

template <class Config, class Greensfunction>
void SpinSusceptibility_X<Config, Greensfunction>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    contrib_dry(UP, UP, func);
    contrib_dry(UP, DOWN, func);
    contrib_dry(DOWN, UP, func);
    contrib_dry(DOWN, DOWN, func);
    return;
}

template <class Config, class Greensfunction>
void SpinSusceptibility_X<Config, Greensfunction>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    GOmegaRetType retval = contrib(UP, UP, configuration, dowick) + contrib(UP, DOWN, configuration, dowick) + contrib(DOWN, UP, configuration, dowick) + contrib(DOWN, DOWN, configuration, dowick);
    this->add_bin(retval/ 4.0);
    return;
}

/**
The Spin-Susceptiblity in Z direction
*/
template <class Config, class Greensfunction>
class SpinSusceptibility_Z : public Network_Cache<Config, typename Greensfunction::GOmegaRetType>
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::Comm Net;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<GFRetVal> Function;
    typedef typename Greensfunction::GOmegaRetType GOmegaRetType;
    typedef GOmegaRetType ObservableType;
    /**
    */
    inline SpinSusceptibility_Z(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "SpinSusceptibility_Z"), len(params.N), beta(params.beta)
    {
        myf = new GOmegaRetType*[8];
        delta_tau = beta/points;
        const uint tablesize = 2 * points;
        for (uint k = 0; k < 8; ++k)
        {
            myf[k] = new GOmegaRetType[tablesize];
            memset(myf[k], 0, tablesize * sizeof(GOmegaRetType));
        }
        for (uint n = 0; n < points; ++n)
        {
	    FPType omegan = M_PI/beta*(2*n+1);
            for (int s = 0; s < 8; ++s)
            {
                SPINS sigma = (s>>2&1? DOWN: UP);
                SPINS sigma_r = (s>>1&1? DOWN: UP);
                SPINS sigma_s = (s&1? DOWN: UP);

                typename Greensfunction::Vertex v1;
                typename Greensfunction::Vertex v2;
                v1.spin = sigma_r;
                v2.spin = sigma;
                GOmegaRetType goma = Greensfunction::gomega(omegan, v1, v2);
                GOmegaRetType gomam = Greensfunction::gomega(-omegan, v1, v2);
                v1.spin = sigma;
                v2.spin = sigma_s;
                GOmegaRetType gomb = Greensfunction::gomega(omegan, v1, v2);
                GOmegaRetType gombm = Greensfunction::gomega(-omegan, v1, v2);
                for (uint k = 0; k < points; ++k)
                {
                    FPType tau = k*delta_tau;
                    std::complex<FPType> expt = exp(std::complex<FPType>(0.0, omegan * tau));
                    std::complex<FPType> expmt = std::conj(expt);
                    accessphi(sigma, sigma_r, sigma_s, k) += goma * gomb * expmt + gomam*gombm*expt;
                    accessphi(sigma, sigma_r, sigma_s, points + k) += goma * gomb * expt + gomam*gombm*expmt;
                }
            }
        }
        for (uint s = 0; s < 8; ++s)
            for (uint k = 0; k < tablesize; ++k)
                myf[s][k] /= beta;
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    This determines the Spin-Susceptiblity for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (GCC_VERSION > GCC_VER(4,5,0))
    static constexpr uint points = 1000;
#else
    static const uint points = 1000;
#endif
    const uint32_t& len;
    FPType beta;
    FPType delta_tau;
    GOmegaRetType** myf;
    GOmegaRetType& accessphi(SPINS sigma, SPINS sigma_r, SPINS sigma_s, int tau_idx)
    {
        return myf[
                   (sigma == UP ? 0: 4) +
                   (sigma_r == UP ? 0: 2) +
                   (sigma_s == UP ? 0 : 1)
               ][tau_idx];
    }
    GOmegaRetType contrib(const SPINS sigma, const SPINS sigmaprime, const Configuration& config, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
    {
        typename Greensfunction::Vertex v1;
        v1.spin = sigmaprime;
        v1.tau = 0;
        typename Greensfunction::Vertex v2;
        v2.spin = sigmaprime;
        v2.tau = 0;
        GOmegaRetType sum1 = dowick(v1, v2);
        v1.spin = sigma;
        v2.spin = sigma;
        GFRetVal gf1 = beta * Greensfunction::eval(v1, v2);
        GOmegaRetType suma = 0.0;
        for (uint r = 0; r < config.size(); ++r)
            for (uint s = 0; s < config.size(); ++s)
            {
                FPType taurs = config[r].tau - config[s].tau;
                int tau_idx = 0;
                if (taurs < 0)
                {
                    tau_idx += points;
                    taurs = -taurs;
                }
                tau_idx += static_cast<int>(trunc(taurs/delta_tau));
                suma += config.matcont(r, s, UP, UP) * accessphi(sigma, UP, UP, tau_idx);
                suma += config.matcont(r, s, DOWN, UP) * accessphi(sigma, DOWN, UP, tau_idx);
                suma += config.matcont(r, s, UP, DOWN) * accessphi(sigma, UP, DOWN, tau_idx);
                suma += config.matcont(r, s, DOWN, DOWN) * accessphi(sigma, DOWN, DOWN, tau_idx);
            }
        sum1 *= (gf1 - suma);
        GOmegaRetType sum2 = accessphi(sigma, sigmaprime, sigmaprime, 0);
        typename Configuration::MatConf::MatType psi(1, static_cast<unsigned int>(2 * config.size()));
        for (uint k = 0; k < config.size(); ++k)
        {
            v1.tau = config[k].tau;
            v2.tau = 0;
            v1.spin = UP;
            v2.spin = sigmaprime;
            psi(0,2*k) = Greensfunction::eval(v1, v2);
            v1.spin = DOWN;
            psi(0, 2*k + 1) = Greensfunction::eval(v1, v2);
        }
        config.matcont.multiplyVectorbyConfiguration_right(psi);
//        psi = psi * config.matcont.mat;
        GOmegaRetType sum3 = 0.0;
        for (uint k = 0; k < config.size(); ++k)
        {
            FPType tau = config[k].tau;
            int tau_idx = points + static_cast<int>(trunc(tau/delta_tau));/*see formula. that way we take the negative sign into consideration*/    
            sum3 += psi(0, 2*k) * accessphi(sigma, UP, sigmaprime, tau_idx);
            sum3 += psi(0, 2*k + 1) * accessphi(sigma, DOWN, sigmaprime, tau_idx);
        }
        GOmegaRetType sum4 = 0.0;
        typename Configuration::MatConf::MatType chi(static_cast<unsigned int>(2 * config.size()), 1);
        for (uint k = 0; k < config.size(); ++k)
        {
            v1.tau = 0;
            v2.tau = config[k].tau;
            v2.spin = UP;
            v1.spin = sigmaprime;
            chi(2*k, 0) = Greensfunction::eval(v1, v2);
            v2.spin = DOWN;
            chi(2*k + 1, 0) = Greensfunction::eval(v1, v2);
        }
        config.matcont.multiplyVectorbyConfiguration_left(chi);
//        chi = config.matcont.mat * chi;
        for (uint k = 0; k < static_cast<uint>(config.size()); ++k)
        {
            FPType tau = config[k].tau;
            int tau_idx = static_cast<int>(trunc(tau/delta_tau));//tau should be positive being straight from a vertex
            sum4 += chi(2*k, 0) * accessphi(sigma, sigmaprime, UP, tau_idx);
            sum4 += chi(2*k + 1, 0) * accessphi(sigma, sigmaprime, DOWN, tau_idx);
        }
        GOmegaRetType sum5 = 0;//we can reuse psi and chi
        for (uint r = 0; r < static_cast<uint>(config.size()); ++r)
            for (uint s = 0; s < static_cast<uint>(config.size()); ++s)
            {
                FPType taurs = config[s].tau - config[r].tau;
                int tau_idx = 0;
                if (taurs < 0)
                {
                    tau_idx += points;
                    taurs = -taurs;
                }
                tau_idx += static_cast<int>(trunc(taurs/delta_tau));
                sum5 += psi(0, 2*r) * chi(2*s, 0)*accessphi(sigma, UP, UP, tau_idx);
                sum5 += psi(0, 2*r+1) * chi(2*s, 0)*accessphi(sigma, DOWN, UP, tau_idx);
                sum5 += psi(0, 2*r) * chi(2*s+1, 0)*accessphi(sigma, UP, DOWN, tau_idx);
                sum5 += psi(0, 2*r+1) * chi(2*s+1, 0)*accessphi(sigma, DOWN, DOWN, tau_idx);
            }
       return sum1 - sum2 - sum3 - sum4 - sum5;//although Davids thesis states this formula differently careful testing reveals this choice of signs(with sum 3 and sum4 negative) to be the right one. Note that the signs of David's thesis can be reestablished by giving psi and chi an additional sign, or by exchanging the order in which the vertices are evaluated. This shouldn't change the sign of sum5 since there's a product of psi and chi and hence any sign cancels.
    }
    void contrib_dry(SPINS sigma, SPINS sigmaprime, DryRun<typename Configuration::value_type, GFRetVal>& func)
    {//yes, the only thing that we don't derive from tables only depends on sigmaprime 
        typename Greensfunction::Vertex v1;
        v1.spin = sigmaprime;
        v1.tau = 0;
        func(v1, v1);
    }
};

template <class Config, class Greensfunction>
void SpinSusceptibility_Z<Config, Greensfunction>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    contrib_dry(UP, UP, func);
    contrib_dry(UP, DOWN, func);
    contrib_dry(DOWN, UP, func);
    contrib_dry(DOWN, DOWN, func);
    return;
}

template <class Config, class Greensfunction>
void SpinSusceptibility_Z<Config, Greensfunction>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    GOmegaRetType retval = contrib(UP, UP, configuration, dowick) - contrib(UP, DOWN, configuration, dowick) - contrib(DOWN, UP, configuration, dowick) + contrib(DOWN, DOWN, configuration, dowick);
    this->add_bin(retval/ 4.0);
    return;
}

template <typename FPType>
struct Omegadata
{
  FPType omega;
  FPType omegasq;
  FPType invomega;
};

/**
A class for measuring the local Greensfunctions on the impurity surrounding bath sites
*/
template <class Config, class GreensFunction, SPINS Spin>
class LocalBathGreensfunctions : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef GreensFunction GF;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///<we pack all 80000 functions into a single valarray...
    /**
    The Constructor for the LocalBathGreensfunctions.
    Notice that in unmixed[]  we tabulate <d^\dagger c>
    */
    LocalBathGreensfunctions(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "LocalBathGreensfunctions"), functionpoints(params.functionpoints), delta_s(params.delta_s),
            beta(params.beta),
            Ny(params.Nb*2),
            Nx(params.Nx),
            Nt(300 + 1),
            delta_beta(params.beta/static_cast<FPType>(Nt - 1)),
            Nw(4*Nt),
            Nyt(Nt * Ny),
            Nxyt(Nyt*params.Nx),
            functions(Ny*params.Nx)
    {
#ifdef _OPENMP
        double start2 = omp_get_wtime();
#endif
        std::cout<<"Creating Local Bath Greensfunction"<<std::endl;
        const GOmegaData<FPType> *const updata = GF::gomegaup->data;
        const GOmegaData<FPType> *const downdata = GF::gomegadown->data;
        //let's set up the Q(r,n,omega_n) data
        Omegadata<FPType>* omega = new Omegadata<FPType>[Nw/2];
	std::complex<FPType>* gup = new std::complex<FPType>[Nw];
	std::complex<FPType>* gdown = new std::complex<FPType>[Nw];
        for (int t = 0; t < Nw/2; ++t)
	{
	    gup[2*t] = conj((*GF::gomegaup)(t));
	    gup[2*t + 1] = conj((*GF::gomegaup)( -t - 1));
	    gdown[2*t] = conj((*GF::gomegadown)(t));
	    gdown[2*t + 1] = conj((*GF::gomegadown)( -t - 1));
            omega[t].omega = M_PI/params.beta*(2*t + 1);
	    omega[t].omegasq = omega[t].omega*omega[t].omega;
	    omega[t].invomega = 1.0/omega[t].omega;
	}
        FPType pref = GF::gethybridization() / std::sqrt(params.Nx);
        unmixeddata = new std::complex<float>[2*Nxyt];
	//Since we have a simple, analytical expression for the tau dependence of the offset of the unmixed greensfunctions
	//we write that one out first to the unmixeddata array
	FPType* cosharrayup = new FPType[functions];
	FPType* cosharraydown = new FPType[functions];
	for(uint k = 0; k < params.Nx; ++k)
	  for(uint m = 0; m < Ny; ++m)
	  {
	    cosharrayup[k*Ny + m]   = 1.0/std::cosh(params.beta*updata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	    cosharraydown[k*Ny + m] = 1.0/std::cosh(params.beta*downdata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	  }
#pragma omp parallel for
	for(uint n = 0; n < Ny; ++n)
	{
	  FPType* tempup = new FPType[functions];
	  FPType* tempdown = new FPType[functions];
	  
	  for(uint k = 0; k < params.Nx; ++k)
	    for(uint m = 0; m < Ny; ++m)
	    {
	      tempup[k*Ny + m] = cosharrayup[k*Ny + m] * norm(updata[k*Ny + /*n*/m].evec[/*m*/n]);
	      tempdown[k*Ny + m] = cosharraydown[k*Ny + m] * norm(downdata[k*Ny + /*n*/m].evec[/*m*/n]);
	    }
	  for(uint t = 0; t < Nt; ++t)
	  {
	    FPType ftup = 0.0;
	    FPType ftdown = 0.0;
	    FPType arg = t * delta_beta - params.beta/2.0;
	    FPType cup = 0.0;
	    FPType cdown = 0.0;
	    for(uint k = 0; k < params.Nx; ++k)//being careful we employ Kahan summation
	      for(uint m = 0; m < Ny; ++m)
	      {
		FPType argup = tempup[k*Ny + m] * std::exp(updata[k*Ny + m].lambda*arg);
		FPType argdown = tempdown[k*Ny + m] * std::exp(downdata[k*Ny + m].lambda*arg);
		FPType y = argup - cup;
		FPType t = ftup + y;
		cup = (t - ftup) - y;
		ftup = t;
		
		y = argdown - cdown;
		t = ftdown + y;
		cdown = (t - ftdown) - y;
		ftdown = t;
	      }
	   for(uint r = 0; r < params.Nx; ++r)
	   {
	      (unmixeddata        + n*Nt*params.Nx + r*Nt)[t] = static_cast<float>(ftup);
              (unmixeddata + Nxyt + n*Nt*params.Nx + r*Nt)[t] = static_cast<float>(ftdown);
	   }
	  }
	  delete [] tempup;
	  delete [] tempdown;
	}
	delete [] cosharrayup;
	delete [] cosharraydown;
        const uint Nxw = Nw*params.Nx;
	const std::complex<FPType> expNt = std::exp(std::complex<FPType>(0.0, -M_PI/(Nt-1)));
	const FPType normierungunmixed = 1.0/params.beta / params.Nx;
	const FPType normierungmixed = 1.0/params.beta / std::sqrt(params.Nx);
	mixeddata = new std::complex<float>[2*Nxyt];
	std::ofstream gimag("gimag.txt");
//#pragma omp parallel for
        for (uint n = 0; n < Ny; ++n)
        {
	    std::complex<FPType>* Qup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* Qdown = new std::complex<FPType>[Nxw];
	    std::complex<FPType>* funup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* fundown = new std::complex<FPType>[Nxw];
            memset(funup, 0, Nxw*sizeof(std::complex<FPType>));
            memset(fundown, 0, Nxw*sizeof(std::complex<FPType>));
//            double start = omp_get_wtime();
            for (uint k = 0; k < params.Nx; ++k)
            {
                for (uint m = 0; m < Ny; ++m)
                {
                    std::complex<FPType> facup   = conj(updata  [k*Ny + m].u) * updata[k*Ny + m].evec[n];
                    std::complex<FPType> facdown = conj(downdata[k*Ny + m].u) * downdata[k*Ny + m].evec[n];
                    for (int omega_idx = 0; omega_idx < Nw/2; ++omega_idx)
                    {//the layout of the frequencies is now (w_n, -w_n) , that is every negative frequency is stored next to its positive counterpart. 
                     //Hopefully this gives a better data locality
                        funup[2*omega_idx*params.Nx + k]         += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        funup[(2*omega_idx + 1)*params.Nx + k]   += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, omega[omega_idx].omega);
                        fundown[2*omega_idx*params.Nx + k]       += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        fundown[(2*omega_idx + 1)*params.Nx + k] += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, omega[omega_idx].omega);
                    }
                }
            }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint w = 0; w < Nw; ++w)
            {
                fourier1(reinterpret_cast<FPType*>(funup + w*params.Nx), params.Nx, 1);
                fourier1(reinterpret_cast<FPType*>(fundown + w*params.Nx), params.Nx, 1);
		//funup as well as fundown now contain Q(r, i omega) for a particular value of the orbital n
                for (uint r = 0; r < params.Nx; ++r)
                {
		    funup[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
		    fundown[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
                    (Qup + r*Nw)[w] = funup[w*params.Nx + r];//norm(funup[w*params.Nx + r]);
                    (Qdown + r*Nw)[w] = fundown[w*params.Nx + r];//norm(fundown[w*params.Nx + r]);
//		    funup[w*params.Nx + r] = funup[w*params.Nx + r];
//		    fundown[w*params.Nx + r] = fundown[w*params.Nx + r];
                }
            }
/*            for(uint r = 0; r < params.Nx; ++r)
	    {
	      for(uint w = 0; w < Nw/2; ++w)
	      {
		gimag<<omega[w].omega<<" "<<real(funup[2*w*params.Nx + r]*gup[2*w])<<std::endl;
	      }
	      gimag<<"&"<<std::endl;
	    }*/
            for(uint r = 0; r < params.Nx; ++r)
            for(uint w = 0; w < Nw/2; ++w)
	    {
	      std::complex<FPType> temp = conj((Qup + r*Nw)[2*w]);
	      (Qup + r*Nw)[2*w] *= conj((Qup + r*Nw)[2*w+1]);
	      (Qup + r*Nw)[2*w+1] *= temp;
	      temp = conj((Qdown + r*Nw)[2*w]);
	      (Qdown + r*Nw)[2*w] *= conj((Qdown + r*Nw)[2*w+1]);
	      (Qdown + r*Nw)[2*w+1] *= temp;
	    }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint r = 0; r < params.Nx; ++r)
	    {
	        std::complex<FPType> expt = 1;
                for (uint t = 0; t < Nt; ++t)//here is the final Matsubara transform
                {
                    std::complex<FPType> tempup = 0;
                    std::complex<FPType> tempdown = 0;
                    std::complex<FPType> tempupmixed = 0;
                    std::complex<FPType> tempdownmixed = 0;
                    FPType tau = t*delta_beta;
		    std::complex<FPType> expiom = expt;//std::exp(tau * std::complex<FPType>(0.0, omega[0].omega));
		    std::complex<FPType> expfac = expiom * expiom;
                    for (int omega_idx = 0; omega_idx < Nw/2; ++omega_idx)
                    {
                        std::complex<FPType> gupp = gup[2*omega_idx];
			std::complex<FPType> gupm = gup[2*omega_idx + 1];
                        std::complex<FPType> gdownp = gdown[2*omega_idx];
			std::complex<FPType> gdownm = gdown[2*omega_idx + 1];
			std::complex<FPType> cexpiom = conj(expiom);
                        tempup += cexpiom * (Qup + r * Nw)[2*omega_idx] * gupp
                               + /*c*/expiom * (Qup + r * Nw)[2*omega_idx + 1] * gupm;
                        tempdown += cexpiom * (Qdown + r * Nw)[2*omega_idx] * gdownp
                                 +/*c*/expiom*(Qdown + r * Nw)[2*omega_idx + 1] * gdownm;

			tempupmixed += cexpiom*(funup[(2*omega_idx)*params.Nx + r] * gupp)
			            + /*c*/expiom * (funup[(2*omega_idx + 1)*params.Nx + r] * gupm);
			tempdownmixed += cexpiom*(fundown[(2*omega_idx)*params.Nx + r] * gdownp)
			              + /*c*/expiom * (fundown[(2*omega_idx + 1)*params.Nx + r] * gdownm);
			expiom *= expfac;
                    }
//                    test<<tempupmixed*normierungmixed<<std::endl;
//The unmixeddata has been debugged in comparison to a straightforward calculation from the real-space Hamiltonian
                    (unmixeddata        + n*Nt*params.Nx + r*Nt)[t] += tempup*normierungunmixed;
                    (unmixeddata + Nxyt + n*Nt*params.Nx + r*Nt)[t] += tempdown*normierungunmixed;
//the sign of the mixeddata is wrong. But since for now mixeddata is only accessed as some squared quantity it doesn't hurt.
                    (mixeddata        + n*Nt*params.Nx + r*Nt)[t] = tempupmixed*normierungmixed;
                    (mixeddata + Nxyt + n*Nt*params.Nx + r*Nt)[t] = tempdownmixed*normierungmixed;
		    expt *= expNt;
                }
	    }
        delete [] funup;
        delete [] fundown;
	delete [] Qup;
	delete [] Qdown;
        }
//        exit(-1);
#ifdef _OPENMP
        std::cout<<"Initialization took "<<omp_get_wtime() - start2<<" seconds"<<std::endl;
#endif
	delete [] gup;
	delete [] gdown;
	delete [] omega;
        std::cout<<"Local Bath GreensFunction done"<<std::endl;
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&) {}
    /**
    this determines the LocalBathGreensfunctions for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& functionpoints;
    const FPType delta_s;
    const FPType beta;
    std::complex<float>* mixeddata;///< here we store <d^+ gamma>. We use floats to keep the memory footprint small
    std::complex<float>* unmixeddata;///< here we store <gamma^+ gamma>
    const uint Ny;
    const uint Nx;
    const uint Nt;
    const FPType delta_beta;
    const uint Nw;
    const uint Nyt;
    const uint Nxyt;
    const uint32_t functions;
    std::complex<FPType> accessmixed(uint r, uint n, SPINS spin, FPType tau1, FPType tau2 ) const
    {
      const FPType tiny = std::numeric_limits<FPType>::epsilon();
      FPType delta_tau = tau1 - tau2;
      FPType sign = 1.0;
      std::complex<float>* dataptr = mixeddata;
      if (spin == DOWN) dataptr += Nxyt;
      dataptr = dataptr + n*Nt*Nx + r*Nt;
      if (std::abs(delta_tau) < tiny)
      {
        //return only the particle number
        return std::complex<FPType>(dataptr[0]);
      }
      if(delta_tau < 0)
      {
	sign = -1.0;
	delta_tau += beta;
      }
      FPType fptau_idx0;
      FPType rem = std::modf(delta_tau/delta_beta, &fptau_idx0);//round to the smaller index and determine how far we're of
      std::size_t tau_idx0 = lround(fptau_idx0);
//    std::cout<<"tau_0: "<<tau_idx0<<" "<<g[tau_idx0]<<std::endl;
      return std::complex<FPType>(lerp(float(rem), dataptr[tau_idx0], dataptr[tau_idx0 + 1]))*sign;
    }
    std::complex<FPType> accessunmixed(uint r, uint n, SPINS spin, FPType tau) const
    {
      std::complex<float>* dataptr = unmixeddata;
      if (spin == DOWN) dataptr += Nxyt;
      dataptr = dataptr + n*Nt*Nx + r*Nt;
      FPType fptau_idx0;
      FPType rem = std::modf(tau/delta_beta, &fptau_idx0);//round to the smaller index and determine how far we're of
      std::size_t tau_idx0 = static_cast<std::size_t>(lround(fptau_idx0));
//    std::cout<<"tau_0: "<<tau_idx0<<" "<<g[tau_idx0]<<std::endl;
//if(tau_idx0 == Nt) return dataptr[tau_idx0];
      return std::complex<FPType>(lerp(float(rem), dataptr[tau_idx0], dataptr[tau_idx0 + 1]));
    }
};

template <class Config, class GreensFunction, SPINS Spin>
void LocalBathGreensfunctions<Config, GreensFunction, Spin>::evaluate(const Configuration& config, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(functions);
#pragma omp parallel for    
    for (unsigned int k = 0; k < func.size(); ++k)
    {
        func[k].resize(functionpoints);
	uint n = k / Nx;
	uint r = k % Nx;
        for (unsigned int j = 0; j < functionpoints/*-1*/; ++j)
        {
            const FPType tau = j * delta_s;
//	    uint unmixedidx = trunc(tau/delta_beta);
            GFRetVal t1 = accessunmixed(r, n, Spin, tau);
//            (unmixeddata        + k*Nt + r*Nt)[unmixedidx];//access the up-sector
            GFRetVal t2 = 0;
	    for (uint q = 0; q < config.size(); ++q)
            for (uint s = 0; s < config.size(); ++s)
            {
	      if(Spin == UP)
	      t2 += config.matcont.mat(2*q, 2*s) * accessmixed(r, n, UP, config[q].tau, 0) * conj(accessmixed(r, n, UP, tau, config[s].tau));
	      t2 += config.matcont.mat(2*q, 2*s + 1) * static_cast<FPType>(0.0);//For now disabled
              t2 += config.matcont.mat(2*q + 1, 2*s) * static_cast<FPType>(0.0);//For now disabled
	      if(Spin == DOWN)
              t2 += config.matcont.mat(2*q + 1, 2*s + 1) * accessmixed(r, n, DOWN, config[q].tau, 0) * conj(accessmixed(r, n, DOWN, tau, config[s].tau));
            }
            
            //add to measurement
            func[k][j] = (t1 - t2)* config.phase;
        }
/*            const FPType tau = beta-0.001;
//	    uint unmixedidx = trunc(tau/delta_beta);
            GFRetVal t1 = accessunmixed(r, n, UP, tau);
//            (unmixeddata        + k*Nt + r*Nt)[unmixedidx];//access the up-sector
            GFRetVal t2 = 0;
	    for (uint r = 0; r < config.size(); ++r)
            for (uint s = 0; s < config.size(); ++s)
            {
	      t2 += config.matcont.mat(2*r, 2*s) * accessmixed(r, n, UP, config[r].tau, 0) * conj(accessmixed(r, n, UP, tau, config[s].tau));
	      t2 += config.matcont.mat(2*r, 2*s + 1) * static_cast<FPType>(0.0);//For now disabled
              t2 += config.matcont.mat(2*r + 1, 2*s) * static_cast<FPType>(0.0);//For now disabled
              t2 += config.matcont.mat(2*r + 1, 2*s + 1) * accessmixed(r, n, DOWN, config[r].tau, 0) * conj(accessmixed(r, n, DOWN, tau, config[s].tau));
            }
            
            //add to measurement
            func[k][functionpoints - 1] = (t1 - t2)* config.phase;*/
    }
    this->add_bin(func);
    return;
}

/**
A class for measuring the local Greensfunctions on the impurity surrounding bath sites, averaged over both spin sectors
*/
template <class Config, class GreensFunction>
class LocalBathGreensfunctions_averaged : public Network_Cache<Config, std::valarray<std::valarray<typename Config::SignType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef GreensFunction GF;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<Function> ObservableType;///<we pack all 80000 functions into a single valarray...
    /**
    The Constructor for the LocalBathGreensfunctions.
    Notice that in unmixed[]  we tabulate <d^\dagger c>
    */
    LocalBathGreensfunctions_averaged(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config, ObservableType>(n, "LocalBathGreensfunctions_averaged"), functionpoints(params.functionpoints), delta_s(params.delta_s),
            beta(params.beta),
            Ny(params.Nb*2),
            Nx(params.Nx),
            Nt(300 + 1),
            delta_beta(params.beta/static_cast<FPType>(Nt - 1)),
            Nw(4*Nt),
            Nyt(Nt * Ny),
            Nxyt(Nyt*params.Nx),
            functions(Ny*params.Nx)
    {
#ifdef _OPENMP
        double start2 = omp_get_wtime();
#endif
        std::cout<<"Creating Local Bath Greensfunction(the averaged one...)"<<std::endl;
        const GOmegaData<FPType> *const updata = GF::gomegaup->data;
        const GOmegaData<FPType> *const downdata = GF::gomegadown->data;
        //let's set up the Q(r,n,omega_n) data
        Omegadata<FPType>* omega = new Omegadata<FPType>[Nw/2];
	std::complex<FPType>* gup = new std::complex<FPType>[Nw];
	std::complex<FPType>* gdown = new std::complex<FPType>[Nw];
        for (int t = 0; t < Nw/2; ++t)
	{
	    gup[2*t] = conj((*GF::gomegaup)(t));
	    gup[2*t + 1] = conj((*GF::gomegaup)( -t - 1));
	    gdown[2*t] = conj((*GF::gomegadown)(t));
	    gdown[2*t + 1] = conj((*GF::gomegadown)( -t - 1));
            omega[t].omega = M_PI/params.beta*(2*t + 1);
	    omega[t].omegasq = omega[t].omega*omega[t].omega;
	    omega[t].invomega = 1.0/omega[t].omega;
	}
        FPType pref = GF::gethybridization() / std::sqrt(params.Nx);
        unmixeddata = new std::complex<float>[2*Nxyt];
	//Since we have a simple, analytical expression for the tau dependence of the offset of the unmixed greensfunctions
	//we write that one out first to the unmixeddata array
	FPType* cosharrayup = new FPType[functions];
	FPType* cosharraydown = new FPType[functions];
	for(uint k = 0; k < params.Nx; ++k)
	  for(uint m = 0; m < Ny; ++m)
	  {
	    cosharrayup[k*Ny + m]   = 1.0/std::cosh(params.beta*updata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	    cosharraydown[k*Ny + m] = 1.0/std::cosh(params.beta*downdata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	  }
#pragma omp parallel for
	for(uint n = 0; n < Ny; ++n)
	{
	  FPType* tempup = new FPType[functions];
	  FPType* tempdown = new FPType[functions];
	  
	  for(uint k = 0; k < params.Nx; ++k)
	    for(uint m = 0; m < Ny; ++m)
	    {
	      tempup[k*Ny + m] = cosharrayup[k*Ny + m] * norm(updata[k*Ny + /*n*/m].evec[/*m*/n]);
	      tempdown[k*Ny + m] = cosharraydown[k*Ny + m] * norm(downdata[k*Ny + /*n*/m].evec[/*m*/n]);
	    }
	  for(uint t = 0; t < Nt; ++t)
	  {
	    FPType ftup = 0.0;
	    FPType ftdown = 0.0;
	    FPType arg = t * delta_beta - params.beta/2.0;
	    FPType cup = 0.0;
	    FPType cdown = 0.0;
	    for(uint k = 0; k < params.Nx; ++k)//being careful we employ Kahan summation
	      for(uint m = 0; m < Ny; ++m)
	      {
		FPType argup = tempup[k*Ny + m] * std::exp(updata[k*Ny + m].lambda*arg);
		FPType argdown = tempdown[k*Ny + m] * std::exp(downdata[k*Ny + m].lambda*arg);
		FPType y = argup - cup;
		FPType t = ftup + y;
		cup = (t - ftup) - y;
		ftup = t;
		
		y = argdown - cdown;
		t = ftdown + y;
		cdown = (t - ftdown) - y;
		ftdown = t;
	      }
	   for(uint r = 0; r < params.Nx; ++r)
	   {
	      (unmixeddata        + n*Nt*params.Nx + r*Nt)[t] = static_cast<float>(ftup);
              (unmixeddata + Nxyt + n*Nt*params.Nx + r*Nt)[t] = static_cast<float>(ftdown);
	   }
	  }
	  delete [] tempup;
	  delete [] tempdown;
	}
	delete [] cosharrayup;
	delete [] cosharraydown;
        const uint Nxw = Nw*params.Nx;
	const std::complex<FPType> expNt = std::exp(std::complex<FPType>(0.0, -M_PI/(Nt-1)));
	const FPType normierungunmixed = 1.0/params.beta / params.Nx;
	const FPType normierungmixed = 1.0/params.beta / std::sqrt(params.Nx);
	mixeddata = new std::complex<float>[2*Nxyt];
#pragma omp parallel for
        for (uint n = 0; n < Ny; ++n)
        {
	    std::complex<FPType>* Qup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* Qdown = new std::complex<FPType>[Nxw];
	    std::complex<FPType>* funup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* fundown = new std::complex<FPType>[Nxw];
            memset(funup, 0, Nxw*sizeof(std::complex<FPType>));
            memset(fundown, 0, Nxw*sizeof(std::complex<FPType>));
//            double start = omp_get_wtime();
            for (uint k = 0; k < params.Nx; ++k)
            {
                for (uint m = 0; m < Ny; ++m)
                {
                    std::complex<FPType> facup   = conj(updata  [k*Ny + m].u) * updata[k*Ny + m].evec[n];
                    std::complex<FPType> facdown = conj(downdata[k*Ny + m].u) * downdata[k*Ny + m].evec[n];
                    for (int omega_idx = 0; omega_idx < Nw/2; ++omega_idx)
                    {//the layout of the frequencies is now (w_n, -w_n) , that is every negative frequency is stored next to its positive counterpart. 
                     //Hopefully this gives a better data locality
                        funup[2*omega_idx*params.Nx + k]         += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        funup[(2*omega_idx + 1)*params.Nx + k]   += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, omega[omega_idx].omega);
                        fundown[2*omega_idx*params.Nx + k]       += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        fundown[(2*omega_idx + 1)*params.Nx + k] += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, omega[omega_idx].omega);
                    }
                }
            }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint w = 0; w < Nw; ++w)
            {
                fourier1(reinterpret_cast<FPType*>(funup + w*params.Nx), params.Nx, 1);
                fourier1(reinterpret_cast<FPType*>(fundown + w*params.Nx), params.Nx, 1);
		//funup as well as fundown now contain Q(r, i omega) for a particular value of the orbital n
                for (uint r = 0; r < params.Nx; ++r)
                {
		    funup[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
		    fundown[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
                    (Qup + r*Nw)[w] = funup[w*params.Nx + r];//norm(funup[w*params.Nx + r]);
                    (Qdown + r*Nw)[w] = fundown[w*params.Nx + r];//norm(fundown[w*params.Nx + r]);
//		    funup[w*params.Nx + r] = funup[w*params.Nx + r];
//		    fundown[w*params.Nx + r] = fundown[w*params.Nx + r];
                }
            }
            for(uint r = 0; r < params.Nx; ++r)
            for(uint w = 0; w < Nw/2; ++w)
	    {
	      std::complex<FPType> temp = conj((Qup + r*Nw)[2*w]);
	      (Qup + r*Nw)[2*w] *= conj((Qup + r*Nw)[2*w+1]);
	      (Qup + r*Nw)[2*w+1] *= temp;
	      temp = conj((Qdown + r*Nw)[2*w]);
	      (Qdown + r*Nw)[2*w] *= conj((Qdown + r*Nw)[2*w+1]);
	      (Qdown + r*Nw)[2*w+1] *= temp;
	    }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint r = 0; r < params.Nx; ++r)
	    {
	        std::complex<FPType> expt = 1;
                for (uint t = 0; t < Nt; ++t)//here is the final Matsubara transform
                {
                    std::complex<FPType> tempup = 0;
                    std::complex<FPType> tempdown = 0;
                    std::complex<FPType> tempupmixed = 0;
                    std::complex<FPType> tempdownmixed = 0;
                    FPType tau = t*delta_beta;
		    std::complex<FPType> expiom = expt;//std::exp(tau * std::complex<FPType>(0.0, omega[0].omega));
		    std::complex<FPType> expfac = expiom * expiom;
                    for (int omega_idx = 0; omega_idx < Nw/2; ++omega_idx)
                    {
                        std::complex<FPType> gupp = gup[2*omega_idx];
			std::complex<FPType> gupm = gup[2*omega_idx + 1];
                        std::complex<FPType> gdownp = gdown[2*omega_idx];
			std::complex<FPType> gdownm = gdown[2*omega_idx + 1];
			std::complex<FPType> cexpiom = conj(expiom);
                        tempup += cexpiom * (Qup + r * Nw)[2*omega_idx] * gupp
                               + /*c*/expiom * (Qup + r * Nw)[2*omega_idx + 1] * gupm;
                        tempdown += cexpiom * (Qdown + r * Nw)[2*omega_idx] * gdownp
                                 +/*c*/expiom*(Qdown + r * Nw)[2*omega_idx + 1] * gdownm;

			tempupmixed += cexpiom*(funup[(2*omega_idx)*params.Nx + r] * gupp)
			            + /*c*/expiom * (funup[(2*omega_idx + 1)*params.Nx + r] * gupm);
			tempdownmixed += cexpiom*(fundown[(2*omega_idx)*params.Nx + r] * gdownp)
			              + /*c*/expiom * (fundown[(2*omega_idx + 1)*params.Nx + r] * gdownm);
			expiom *= expfac;
                    }
//                    test<<tempupmixed*normierungmixed<<std::endl;
//The unmixeddata has been debugged in comparison to a straightforward calculation from the real-space Hamiltonian
                    (unmixeddata        + n*Nt*params.Nx + r*Nt)[t] += tempup*normierungunmixed;
                    (unmixeddata + Nxyt + n*Nt*params.Nx + r*Nt)[t] += tempdown*normierungunmixed;
//the sign of the mixeddata is wrong. But since for now mixeddata is only accessed as some squared quantity it doesn't hurt.
                    (mixeddata        + n*Nt*params.Nx + r*Nt)[t] = tempupmixed*normierungmixed;
                    (mixeddata + Nxyt + n*Nt*params.Nx + r*Nt)[t] = tempdownmixed*normierungmixed;
		    expt *= expNt;
                }
	    }
        delete [] funup;
        delete [] fundown;
	delete [] Qup;
	delete [] Qdown;
        }
#ifdef _OPENMP
        std::cout<<"Initialization took "<<omp_get_wtime() - start2<<" seconds"<<std::endl;
#endif
	delete [] gup;
	delete [] gdown;
	delete [] omega;
        std::cout<<"Local Bath GreensFunction done"<<std::endl;
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&) {}
    /**
    This determines the LocalBathGreensfunctions for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& functionpoints;
    const FPType delta_s;
    const FPType beta;
    std::complex<float>* mixeddata;///< here we store <d^+ gamma>. We use floats to keep the memory footprint small
    std::complex<float>* unmixeddata;///< here we store <gamma^+ gamma>
    const uint Ny;
    const uint Nx;
    const uint Nt;
    const FPType delta_beta;
    const uint Nw;
    const uint Nyt;
    const uint Nxyt;
    const uint32_t functions;
    std::complex<FPType> accessmixed(uint r, uint n, SPINS spin, FPType tau1, FPType tau2 ) const
    {
      const FPType tiny = std::numeric_limits<FPType>::epsilon();
      FPType delta_tau = tau1 - tau2;
      FPType sign = 1.0;
      std::complex<float>* dataptr = mixeddata;
      if (spin == DOWN) dataptr += Nxyt;
      dataptr = dataptr + n*Nt*Nx + r*Nt;
      if (std::abs(delta_tau) < tiny)
      {
        //return only the particle number
        return std::complex<FPType>(dataptr[0]);
      }
      if(delta_tau < 0)
      {
	sign = -1.0;
	delta_tau += beta;
      }
      FPType fptau_idx0;
      FPType rem = std::modf(delta_tau/delta_beta, &fptau_idx0);//round to the smaller index and determine how far we're of
      std::size_t tau_idx0 = lround(fptau_idx0);
//    std::cout<<"tau_0: "<<tau_idx0<<" "<<g[tau_idx0]<<std::endl;
      return std::complex<FPType>(lerp(float(rem), dataptr[tau_idx0], dataptr[tau_idx0 + 1]))*sign;
    }
    std::complex<FPType> accessunmixed(uint r, uint n, SPINS spin, FPType tau) const
    {
      std::complex<float>* dataptr = unmixeddata;
      if (spin == DOWN) dataptr += Nxyt;
      dataptr = dataptr + n*Nt*Nx + r*Nt;
      FPType fptau_idx0;
      FPType rem = std::modf(tau/delta_beta, &fptau_idx0);//round to the smaller index and determine how far we're of
      std::size_t tau_idx0 = static_cast<std::size_t>(lround(fptau_idx0));
//    std::cout<<"tau_0: "<<tau_idx0<<" "<<g[tau_idx0]<<std::endl;
//if(tau_idx0 == Nt) return dataptr[tau_idx0];
      return std::complex<FPType>(lerp(float(rem), dataptr[tau_idx0], dataptr[tau_idx0 + 1]));
    }
};

template <class Config, class GreensFunction>
void LocalBathGreensfunctions_averaged<Config, GreensFunction>::evaluate(const Configuration& config, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(functions);
//    ofstream file("t2.txt");
//#pragma omp parallel for    
    for (unsigned int k = 0; k < func.size(); ++k)
    {
        func[k].resize(functionpoints);
	uint n = k / Nx;
	uint r = k % Nx;
        for (unsigned int j = 0; j < functionpoints/*-1*/; ++j)
        {
            const FPType tau = j * delta_s;
//	    uint unmixedidx = trunc(tau/delta_beta);
            GFRetVal t1 = accessunmixed(r, n, UP, tau) + accessunmixed(r, n, DOWN, tau);
//            (unmixeddata        + k*Nt + r*Nt)[unmixedidx];//access the up-sector
            GFRetVal t2 = 0;
	    for (uint q = 0; q < config.size(); ++q)
            for (uint s = 0; s < config.size(); ++s)
            {
	      t2 += config.matcont.mat(2*q, 2*s) * accessmixed(r, n, UP, config[q].tau, 0) * conj(accessmixed(r, n, UP, config[s].tau, tau));
//	      t2 += config.matcont.mat(2*q, 2*s + 1) * static_cast<FPType>(0.0);//For now disabled
//              t2 += config.matcont.mat(2*q + 1, 2*s) * static_cast<FPType>(0.0);//For now disabled
              t2 += config.matcont.mat(2*q + 1, 2*s + 1) * accessmixed(r, n, DOWN, config[q].tau, 0) * conj(accessmixed(r, n, DOWN, config[s].tau, tau));
            }
//            file<<j<<" "<<real(t2)<<std::endl;
            //add to measurement
            func[k][j] = (t1 + t2)* config.phase/2.0;
        }
//        file<<"&"<<std::endl;
    }
//    exit(-1);
    this->add_bin(func);
    return;
}

/**
A class for measuring the Kondo-cloud as evidenced by the correlation function
<S^z_d S^z_c (x)>
where d denotes the dot electron and c the bath electron. All this as a function of distance from the dot.
*/
template <class Config, class GreensFunction, SPINS Spin>
class KondoCloud_Z : public Network_Cache<Config, std::valarray<typename Config::SignType> >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef GreensFunction GF;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<GFRetVal> ObservableType;///<the Kondo Cloud has no time-dependence. it depends only on the position
    /**
    The Constructor for the KondoCloud
    */
    KondoCloud_Z(typename Config::Comm& n, const Parameters& params)/* throw()*/ : Network_Cache<Config, ObservableType>(n, "KondoCloud_Z"), functionpoints(params.functionpoints), beta(params.beta), Ny(params.Nb*2),
            Nx(params.Nx),
            Nt(300 + 1),
            Nw(4*Nt),
            delta_beta(params.beta/static_cast<FPType>(Nt-1)),
            sites(Ny*params.Nx),
	    Nxyt(sites*Nt)
    {
#ifdef _OPENMP
        double start2 = omp_get_wtime();
#endif
        std::cout<<"Creating KondoCloud_Z"<<std::endl;
        const GOmegaData<FPType> *const updata = GF::gomegaup->data;
        const GOmegaData<FPType> *const downdata = GF::gomegadown->data;
        //let's set up the Q(r, n, omega_n) data
        Omegadata<FPType>* omega = new Omegadata<FPType>[Nw/2];
	std::complex<FPType>* gup = new std::complex<FPType>[Nw];
	std::complex<FPType>* gdown = new std::complex<FPType>[Nw];
        for (int t = 0; t < static_cast<int>(Nw/2); ++t)
	{
	    gup[2*t] = conj((*GF::gomegaup)(t));
	    gup[2*t + 1] = conj((*GF::gomegaup)( -t - 1));
	    gdown[2*t] = conj((*GF::gomegadown)(t));
	    gdown[2*t + 1] = conj((*GF::gomegadown)( -t - 1));
            omega[t].omega = M_PI/params.beta*(2*t + 1);
	    omega[t].omegasq = omega[t].omega*omega[t].omega;
	    omega[t].invomega = 1.0/omega[t].omega;
	}
        FPType pref = GF::gethybridization() / std::sqrt(params.Nx);
        unmixeddata = new std::complex<FPType>[2*sites];
	//Since we have a simple, analytical expression for the tau dependence of the offset of the unmixed greensfunctions
	//we write that one out first to the unmixeddata array
	FPType* cosharrayup = new FPType[sites];
	FPType* cosharraydown = new FPType[sites];
	for(uint k = 0; k < params.Nx; ++k)
	  for(uint m = 0; m < Ny; ++m)
	  {
	    cosharrayup[k*Ny + m]   = 1.0/std::cosh(params.beta*updata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	    cosharraydown[k*Ny + m] = 1.0/std::cosh(params.beta*downdata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	  }
#pragma omp parallel for
	for(uint n = 0; n < Ny; ++n)
	{
	  FPType* tempup = new FPType[sites];
	  FPType* tempdown = new FPType[sites];
	  for(uint k = 0; k < params.Nx; ++k)
	    for(uint m = 0; m < Ny; ++m)
	    {
	      tempup[k*Ny + m] = cosharrayup[k*Ny + m] * norm(updata[k*Ny + /*n*/m].evec[/*m*/n]);
	      tempdown[k*Ny + m] = cosharraydown[k*Ny + m] * norm(downdata[k*Ny + /*n*/m].evec[/*m*/n]);
	    }
	    FPType ftup = 0.0;
	    FPType ftdown = 0.0;
	    FPType arg = - params.beta/2.0;//from that we only need a tau == 0 quantity
	    FPType cup = 0.0;
	    FPType cdown = 0.0;
	    for(uint k = 0; k < params.Nx; ++k)//being careful we employ Kahan summation
	      for(uint m = 0; m < Ny; ++m)
	      {
		FPType argup = tempup[k*Ny + m] * std::exp(updata[k*Ny + m].lambda*arg);
		FPType argdown = tempdown[k*Ny + m] * std::exp(downdata[k*Ny + m].lambda*arg);
		FPType y = argup - cup;
		FPType t = ftup + y;
		cup = (t - ftup) - y;
		ftup = t;
		y = argdown - cdown;
		t = ftdown + y;
		cdown = (t - ftdown) - y;
		ftdown = t;
	      }
	   for(uint r = 0; r < params.Nx; ++r)
	   {
	      (unmixeddata        + n*params.Nx)[r] = static_cast<float>(ftup);
              (unmixeddata + sites + n*params.Nx)[r] = static_cast<float>(ftdown);
	   }
	  delete [] tempup;
	  delete [] tempdown;
	}

	delete [] cosharrayup;
	delete [] cosharraydown;
        const unsigned int Nxw = Nw*params.Nx;
	FPType normierungunmixed = 1.0/params.beta / params.Nx;
	FPType normierungmixed = 1.0/params.beta / std::sqrt(params.Nx);
	const std::complex<FPType> expNt = std::exp(std::complex<FPType>(0.0, -M_PI/(Nt-1)));
	mixeddata = new std::complex<FPType>[2*Nxyt];
#pragma omp parallel for
        for (uint n = 0; n < Ny; ++n)
        {
	    std::complex<FPType>* Qup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* Qdown = new std::complex<FPType>[Nxw];
	    std::complex<FPType>* funup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* fundown = new std::complex<FPType>[Nxw];
            memset(funup, 0, Nxw*sizeof(std::complex<FPType>));
            memset(fundown, 0, Nxw*sizeof(std::complex<FPType>));
//            double start = omp_get_wtime();
//            std::cout<<"n = "<<n<<std::endl;
            for (uint k = 0; k < params.Nx; ++k)
            {
                for (uint m = 0; m < Ny; ++m)
                {
                    std::complex<FPType> facup   = conj(updata  [k*Ny + m].u) * updata[k*Ny + m].evec[n];
                    std::complex<FPType> facdown = conj(downdata[k*Ny + m].u) * downdata[k*Ny + m].evec[n];
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {//the layout of the frequencies is now (w_n, -w_n) , that is every negative frequency is stored next to its positive counterpart. 
                     //Hopefully this gives a better data locality
                        funup[2*omega_idx*params.Nx + k]         += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        funup[(2*omega_idx + 1)*params.Nx + k]   += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, omega[omega_idx].omega);
                        fundown[2*omega_idx*params.Nx + k]       += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        fundown[(2*omega_idx + 1)*params.Nx + k] += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, omega[omega_idx].omega);
                    }
                }
            }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint w = 0; w < Nw; ++w)
            {
                fourier1(reinterpret_cast<FPType*>(funup + w*params.Nx), params.Nx, 1);
                fourier1(reinterpret_cast<FPType*>(fundown + w*params.Nx), params.Nx, 1);
		//funup as well as fundown now contain Q(r, i omega) for a particular value of the orbital n
                for (uint r = 0; r < params.Nx; ++r)
                {
		    funup[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
		    fundown[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
                    (Qup + r*Nw)[w] = funup[w*params.Nx + r];//norm(funup[w*params.Nx + r]);
                    (Qdown + r*Nw)[w] = fundown[w*params.Nx + r];//norm(fundown[w*params.Nx + r]);
//		    funup[w*params.Nx + r] = funup[w*params.Nx + r];
//		    fundown[w*params.Nx + r] = fundown[w*params.Nx + r];
                }
            }
            for(uint r = 0; r < params.Nx; ++r)
            for(uint w = 0; w < Nw/2; ++w)
	    {
	      std::complex<FPType> temp = conj((Qup + r*Nw)[2*w]);
	      (Qup + r*Nw)[2*w] *= conj((Qup + r*Nw)[2*w+1]);
	      (Qup + r*Nw)[2*w+1] *= temp;
	      temp = conj((Qdown + r*Nw)[2*w]);
	      (Qdown + r*Nw)[2*w] *= conj((Qdown + r*Nw)[2*w+1]);
	      (Qdown + r*Nw)[2*w+1] *= temp;
	    }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint r = 0; r < params.Nx; ++r)
	    {
	      std::complex<FPType> expt = 1;
	      {//an empty block for the unmixeddata
                    std::complex<FPType> tempup = 0;
                    std::complex<FPType> tempdown = 0;
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {
                        std::complex<FPType> gupp = gup[2*omega_idx];
			std::complex<FPType> gupm = gup[2*omega_idx + 1];
                        std::complex<FPType> gdownp = gdown[2*omega_idx];
			std::complex<FPType> gdownm = gdown[2*omega_idx + 1];
                        tempup += (Qup + r * Nw)[2*omega_idx] * gupp + (Qup + r * Nw)[2*omega_idx + 1] * gupm;
                        tempdown += (Qdown + r * Nw)[2*omega_idx] * gdownp + (Qdown + r * Nw)[2*omega_idx + 1] * gdownm;
                    }
                    (unmixeddata        + n*params.Nx)[r] += tempup*normierungunmixed;
                    (unmixeddata + sites + n*params.Nx)[r] += tempdown*normierungunmixed;
	      }
	      for(uint t = 0; t < Nt; ++t)// here is the final Matsubara transform
	      {
                    std::complex<FPType> tempupmixed = 0;
                    std::complex<FPType> tempdownmixed = 0;
		    std::complex<FPType> expiom = expt;
		    std::complex<FPType> expfac = expiom*expiom;
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {
                        std::complex<FPType> gupp = gup[2*omega_idx];
			std::complex<FPType> gupm = gup[2*omega_idx + 1];
                        std::complex<FPType> gdownp = gdown[2*omega_idx];
			std::complex<FPType> gdownm = gdown[2*omega_idx + 1];
			std::complex<FPType> cexpiom = conj(expiom);
			tempupmixed += cexpiom * funup[(2*omega_idx)*params.Nx + r] * gupp + expiom * funup[(2*omega_idx + 1)*params.Nx + r] * gupm;
			tempdownmixed += cexpiom * fundown[(2*omega_idx)*params.Nx + r] * gdownp + expiom * fundown[(2*omega_idx + 1)*params.Nx + r] * gdownm;
			expiom *= expfac;
                    }
                    (mixeddata        + n*params.Nx*Nt + r*Nt)[t] = tempupmixed*normierungmixed;
                    (mixeddata + Nxyt + n*params.Nx*Nt + r*Nt)[t] = tempdownmixed*normierungmixed;
		    expt *= expNt;
	      }
//                if (r > 3 )exit(-1);
//		test<<"&"<<std::endl;
	    }
//        std::cout<<"Initialization took "<<omp_get_wtime() - start<<" seconds"<<std::endl;
        delete [] funup;
        delete [] fundown;
	delete [] Qup;
	delete [] Qdown;
        }
#ifdef _OPENMP
        std::cout<<"Initialization took "<<omp_get_wtime() - start2<<" seconds"<<std::endl;
#endif
	delete [] gup;
	delete [] gdown;
	delete [] omega;
        std::cout<<"KondoCloud_Z done"<<std::endl;
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func) 
    {
    contrib_dry(UP, UP, func);
    contrib_dry(UP, DOWN, func);
    contrib_dry(DOWN, UP, func);
    contrib_dry(DOWN, DOWN, func);      
    }
    /**
    this determines the KondoCloud for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& functionpoints;
    const FPType beta;
    std::complex<FPType>* mixeddata;///< here we store <d^+ gamma>.
    std::complex<FPType>* unmixeddata;///< here we store <gamma^+ gamma>
    const uint Ny;
    const uint Nx;
    const uint Nt;
    const uint Nw;
    const FPType delta_beta;
    const uint32_t sites;
    const uint Nxyt;
    std::complex<FPType> accessmixed(uint r, uint n, SPINS spin1, FPType tau1 , SPINS spin2, FPType tau2) const
    {
      FPType delta_tau = tau1 - tau2;
      FPType sign = 1.0;
      std::complex<FPType>* dataptr = mixeddata;
      if(spin1 != spin2) return 0.0;//FIXME!!!!!!!!!!!!!!!!!!!! only the case for the spin symmetric case!
      if (spin1 == DOWN) dataptr += Nxyt;
      dataptr = dataptr + n*Nx*Nt + r*Nt;
      //if(std::abs(delta_tau) < std::numeric_limits<FPType>::epsilon())
      if(fpequal(tau1, tau2))
	return std::complex<FPType>(dataptr[0]);
      if(delta_tau < 0)
      {
	sign = -1.0;
	delta_tau += beta;
      }
      FPType fptau_idx0;
      FPType rem = std::modf(delta_tau/delta_beta, &fptau_idx0);//round to the smaller index and determine how far we're off.
      std::size_t tau_idx0 = lround(fptau_idx0);
      return std::complex<FPType>(lerp(rem, dataptr[tau_idx0], dataptr[tau_idx0 + 1]))*sign;
    }
    std::complex<FPType> accessunmixed(uint r, uint n, SPINS spin) const
    {
      std::complex<FPType>* dataptr = unmixeddata;
      if (spin == DOWN) dataptr += sites;
      dataptr = dataptr + n*Nx + r;
      return dataptr[0];
    }
    void contrib_dry(SPINS sigma, SPINS sigmaprime, DryRun<typename Configuration::value_type, GFRetVal>& func) const
    {//yes, the only thing that we don't derive from tables only depends on sigmaprime 
        typename GreensFunction::Vertex v1;
        v1.spin = sigmaprime;
        v1.tau = 0;
        func(v1, v1);
    }
    GFRetVal dotghelper(FPType tau1, SPINS spin1, FPType tau2, SPINS spin2, const Configuration& config)
    {
//    auto gdot = [](FPType tau1, SPINS spin1, FPType tau2, SPINS spin2){return GreensFunction::eval(typename GreensFunction::Vertex(tau1, spin1), typename GreensFunction::Vertex(tau2, spin2));};
    struct GDot{GFRetVal operator()
    (FPType tau1, SPINS spin1, FPType tau2, SPINS spin2){return GreensFunction::eval(typename GreensFunction::Vertex(tau1, spin1), typename GreensFunction::Vertex(tau2, spin2));}
    } gdot;
    GFRetVal retval = gdot(tau1, spin1, tau2, spin2);
    for(uint r = 0; r < config.size(); ++r)
    {
      FPType taur = config[r].tau;
      GFRetVal gtaur_tau1_up = gdot(taur, UP, tau1, spin1);
      GFRetVal gtaur_tau1_down = gdot(taur, DOWN, tau1, spin1);
      for(uint s = 0; s < config.size(); ++s)
      {
	FPType taus = config[s].tau;
	GFRetVal gtau2_taus_up = gdot(tau2, spin2, taus, UP);
	GFRetVal gtau2_taus_down = gdot(tau2, spin2, taus, DOWN);
	retval -= gtaur_tau1_up * config.matcont.mat(2*r, 2*s) * gtau2_taus_up;
	retval -= gtaur_tau1_down * config.matcont.mat(2*r+1, 2*s) * gtau2_taus_up;
	retval -= gtaur_tau1_up * config.matcont.mat(2*r, 2*s+1) * gtau2_taus_down;
	retval -= gtaur_tau1_down * config.matcont.mat(2*r+1, 2*s+1) * gtau2_taus_down;
      }
    }
      return retval;
    }
};

template <class Config, class GreensFunction, SPINS Spin>
void KondoCloud_Z<Config, GreensFunction, Spin>::evaluate(const Configuration& config, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(sites);
    typename GreensFunction::Vertex v1(0.0, UP);
    typename GreensFunction::Vertex v2(0.0, DOWN);
    GFRetVal dotcontrib = dowick(v1, v1) - dowick(v2,v2);
    GFRetVal *const dotgf = new GFRetVal[config.size()*2*2*2];
    for(uint i = 0; i < config.size(); ++i)//here we sum up < d^+ d(t) >
    {
      FPType taui = config[i].tau;
	dotgf[4*i + 0] = dotghelper(0, UP, taui, UP, config);
	dotgf[4*i + 1] = dotghelper(0, UP, taui, DOWN, config);
	dotgf[4*i + 2] = dotghelper(0, DOWN, taui, UP, config);
	dotgf[4*i + 3] = dotghelper(0, DOWN, taui, DOWN, config);

	dotgf[4*(i+config.size()) + 0] = dotghelper(taui, UP, 0, UP, config);
	dotgf[4*(i+config.size()) + 1] = dotghelper(taui, UP, 0, DOWN, config);
	dotgf[4*(i+config.size()) + 2] = dotghelper(taui, DOWN, 0, UP, config);
	dotgf[4*(i+config.size()) + 3] = dotghelper(taui, DOWN, 0, DOWN, config);
    }
//    std::ofstream kc("kc.txt");
#pragma omp parallel for
    for (unsigned int k = 0; k < func.size(); ++k)
    {
	uint n = k / Nx;
	uint r = k % Nx;
	GFRetVal unmixedup = accessunmixed(r, n, UP);
	GFRetVal unmixeddown = accessunmixed(r, n, DOWN);
	GFRetVal mixedupup =  accessmixed(r, n, UP, 0.0, UP, 0.0);
	GFRetVal mixedupdown = accessmixed(r, n, UP, 0.0, DOWN, 0.0);
	GFRetVal mixeddownup = accessmixed(r, n, DOWN, 0.0, UP, 0.0);
	GFRetVal mixeddowndown = accessmixed(r, n, DOWN, 0.0, DOWN, 0.0);

	GFRetVal mixedupupconj =  conj(accessmixed(r, n, UP, 0.0, UP, 0.0));
	GFRetVal mixedupdownconj = conj(accessmixed(r, n, UP, 0.0, DOWN, 0.0));
	GFRetVal mixeddownupconj = conj(accessmixed(r, n, DOWN, 0.0, UP, 0.0));
	GFRetVal mixeddowndownconj = conj(accessmixed(r, n, DOWN, 0.0, DOWN, 0.0));
	for(uint q = 0; q < config.size(); ++q)
	{
	  FPType tauq =  config[q].tau;
	  GFRetVal gmixed_rn_tauq_UP_UP = accessmixed(r, n, UP, tauq, UP, 0.0);
	  GFRetVal gmixed_rn_tauq_DOWN_UP = accessmixed(r, n, DOWN, tauq, UP, 0.0);
	  GFRetVal gmixed_rn_tauq_UP_DOWN = accessmixed(r, n, UP, tauq, DOWN, 0.0);
	  GFRetVal gmixed_rn_tauq_DOWN_DOWN = accessmixed(r, n, DOWN, tauq, DOWN, 0.0);
	  GFRetVal dotq0 = dotgf[4*(q+config.size()) + 0];
	  GFRetVal dotq1 = dotgf[4*(q+config.size()) + 1];
	  GFRetVal dotq2 = dotgf[4*(q+config.size()) + 2];
	  GFRetVal dotq3 = dotgf[4*(q+config.size()) + 3];
	  for(uint s = 0; s < config.size(); ++s)
	  {
	    FPType taus =  config[s].tau;
	    GFRetVal gmixed_rn_taus_UP_UP = conj(accessmixed(r, n, UP, taus, UP, 0.0));
	    GFRetVal gmixed_rn_taus_UP_DOWN = conj(accessmixed(r, n, UP, taus, DOWN, 0.0));
	    GFRetVal gmixed_rn_taus_DOWN_UP = conj(accessmixed(r, n, DOWN, taus, UP, 0.0));
	    GFRetVal gmixed_rn_taus_DOWN_DOWN = conj(accessmixed(r, n, DOWN, taus, DOWN, 0.0));
	    GFRetVal matqs = config.matcont.mat(2*q, 2*s);
	    GFRetVal matqsp = config.matcont.mat(2*q, 2*s + 1);
	    GFRetVal matqps = config.matcont.mat(2*q + 1, 2*s);
	    GFRetVal matqpsp = config.matcont.mat(2*q + 1, 2*s + 1);
	    GFRetVal dots0 = dotgf[4*s+0];
	    GFRetVal dots1 = dotgf[4*s+1];
	    GFRetVal dots2 = dotgf[4*s+2];
	    GFRetVal dots3 = dotgf[4*s+3];
	    unmixedup -= (
	      matqs  *    gmixed_rn_tauq_UP_UP  * gmixed_rn_taus_UP_UP
	     +matqsp*    gmixed_rn_tauq_UP_UP  * gmixed_rn_taus_UP_DOWN
	     +matqps*    gmixed_rn_tauq_DOWN_UP* gmixed_rn_taus_UP_UP
	     +matqpsp  *gmixed_rn_tauq_DOWN_UP* gmixed_rn_taus_UP_DOWN
	    );
	    unmixeddown -= (
	      matqs  *    gmixed_rn_tauq_UP_DOWN  * gmixed_rn_taus_DOWN_UP
	     +matqsp*    gmixed_rn_tauq_UP_DOWN  * gmixed_rn_taus_DOWN_DOWN
	     +matqps*    gmixed_rn_tauq_DOWN_DOWN* gmixed_rn_taus_DOWN_UP
	     +matqpsp  *gmixed_rn_tauq_DOWN_DOWN* gmixed_rn_taus_DOWN_DOWN
	    );
	    
	    mixedupup -= (
	      matqs  *    gmixed_rn_tauq_UP_UP  * dots0
	     +matqsp*    gmixed_rn_tauq_UP_UP  * dots1
	     +matqps*    gmixed_rn_tauq_DOWN_UP* dots0
	     +matqpsp  *gmixed_rn_tauq_DOWN_UP* dots1
	    );

	    mixedupdown -= (
	      matqs  *    gmixed_rn_tauq_UP_DOWN  * dots0
	     +matqsp*    gmixed_rn_tauq_UP_DOWN  * dots1
	     +matqps*    gmixed_rn_tauq_DOWN_DOWN* dots0
	     +matqpsp  *gmixed_rn_tauq_DOWN_DOWN* dots1
	    );
	    
	     mixeddownup -= (
	      matqs  *    gmixed_rn_tauq_UP_UP  * dots2
	     +matqsp*    gmixed_rn_tauq_UP_UP  * dots3
	     +matqps*    gmixed_rn_tauq_DOWN_UP* dots2
	     +matqpsp  *gmixed_rn_tauq_DOWN_UP* dots3
	    );

   	    mixeddowndown -= (
	      matqs  *    gmixed_rn_tauq_UP_DOWN  * dots2
	     +matqsp*    gmixed_rn_tauq_UP_DOWN  * dots3
	     +matqps*    gmixed_rn_tauq_DOWN_DOWN* dots2
	     +matqpsp  *gmixed_rn_tauq_DOWN_DOWN* dots3
	    );
	    
	    GFRetVal gmixed_rn_UP_UP_taus = conj(accessmixed(r, n, UP, 0.0, UP, taus));
	    GFRetVal gmixed_rn_UP_DOWN_taus = conj(accessmixed(r, n, UP, 0.0, DOWN, taus));
	    mixedupupconj -= (
	      matqs   * gmixed_rn_UP_UP_taus  * dotq0
	     +matqsp  * gmixed_rn_UP_DOWN_taus  * dotq0
	     +matqps  * gmixed_rn_UP_UP_taus* dotq2
	     +matqpsp * gmixed_rn_UP_DOWN_taus* dotq2
	    );

	    mixedupdownconj -= (
	      matqs   * gmixed_rn_UP_UP_taus  * dotq1
	     +matqsp  * gmixed_rn_UP_DOWN_taus  * dotq1
	     +matqps  * gmixed_rn_UP_UP_taus* dotq3
	     +matqpsp * gmixed_rn_UP_DOWN_taus* dotq3
	    );
	    
	    GFRetVal gmixed_rn_DOWN_DOWN_taus = conj(accessmixed(r, n, DOWN, 0.0, DOWN, taus));
	    GFRetVal gmixed_rn_DOWN_UP_taus = conj(accessmixed(r, n, DOWN, 0.0, UP, taus));
	    
	     mixeddownupconj -= (
	      matqs   * gmixed_rn_DOWN_UP_taus  * dotq0
	     +matqsp  * gmixed_rn_DOWN_DOWN_taus  * dotq0
	     +matqps  * gmixed_rn_DOWN_UP_taus* dotq2
	     +matqpsp * gmixed_rn_DOWN_DOWN_taus* dotq2
	    );

   	    mixeddowndownconj -= (
	      matqs   * gmixed_rn_DOWN_UP_taus  * dotq1
	     +matqsp  * gmixed_rn_DOWN_DOWN_taus  * dotq1
	     +matqps  * gmixed_rn_DOWN_UP_taus* dotq3
	     +matqpsp * gmixed_rn_DOWN_DOWN_taus* dotq3
	    );
	  }
	}
            //add to measurement
            func[k] = (dotcontrib * (unmixedup - unmixeddown) - mixedupup*mixedupupconj 
            /*+ mixedupdown*mixedupdownconj + mixeddownup*mixeddownupconj*/ 
	    - mixeddowndown*mixeddowndownconj)* config.phase/4.0;
//kc<<dotcontrib <<" "<< (unmixedup - unmixeddown)<<" "<<mixedupup<<" "<<mixedupupconj<<" "<< mixedupdown<<" "<<mixedupdownconj <<" "<< mixeddownup<<" "<<mixeddownupconj <<" "<< mixeddowndown<<" "<<mixeddowndownconj<<std::endl;
    }
//    exit(-1);
    this->add_bin(func);
    delete [] dotgf;
    return;
}

/**
A class for measuring the Kondo-cloud as evidenced by the correlation function
<S^x_d S^x_c (r)>
where d denotes the dot electron and c the bath electron. All this as a function of distance r from the dot.
*/
template <class Config, class GreensFunction, SPINS Spin>
class KondoCloud_X : public Network_Cache<Config, std::valarray<typename Config::SignType> >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef GreensFunction GF;
    typedef std::valarray<GFRetVal> Function;
    typedef std::valarray<GFRetVal> ObservableType;///<the Kondo Cloud has no time-dependence. it depends only on the position
    /**
    The Constructor for the KondoCloud measured along the X- direction
    */
    KondoCloud_X(typename Config::Comm& n, const Parameters& params)/* throw()*/ : Network_Cache<Config, ObservableType>(n, "KondoCloud_X"), functionpoints(params.functionpoints), beta(params.beta), Ny(params.Nb*2),
            Nx(params.Nx),
            Nt(300 + 1),
            Nw(4*Nt),
            delta_beta(params.beta/static_cast<FPType>(Nt-1)),
            sites(Ny*params.Nx),
	    Nxyt(sites*Nt)
    {
#ifdef _OPENMP
        double start2 = omp_get_wtime();
#endif
        std::cout<<"Creating KondoCloud_X"<<std::endl;
        const GOmegaData<FPType> *const updata = GF::gomegaup->data;
        const GOmegaData<FPType> *const downdata = GF::gomegadown->data;
        //let's set up the Q(r, n, omega_n) data
        Omegadata<FPType>* omega = new Omegadata<FPType>[Nw/2];
	std::complex<FPType>* gup = new std::complex<FPType>[Nw];
	std::complex<FPType>* gdown = new std::complex<FPType>[Nw];
        for (int t = 0; t < static_cast<int>(Nw/2); ++t)
	{
	    gup[2*t] = conj((*GF::gomegaup)(t));
	    gup[2*t + 1] = conj((*GF::gomegaup)( -t - 1));
	    gdown[2*t] = conj((*GF::gomegadown)(t));
	    gdown[2*t + 1] = conj((*GF::gomegadown)( -t - 1));
            omega[t].omega = M_PI/params.beta*(2*t + 1);
	    omega[t].omegasq = omega[t].omega*omega[t].omega;
	    omega[t].invomega = 1.0/omega[t].omega;
	}
        FPType pref = GF::gethybridization() / std::sqrt(params.Nx);
        unmixeddata = new std::complex<FPType>[2*sites];
	//Since we have a simple, analytical expression for the tau dependence of the offset of the unmixed greensfunctions
	//we write that one out first to the unmixeddata array
	FPType* cosharrayup = new FPType[sites];
	FPType* cosharraydown = new FPType[sites];
	for(uint k = 0; k < params.Nx; ++k)
	  for(uint m = 0; m < Ny; ++m)
	  {
	    cosharrayup[k*Ny + m]   = 1.0/std::cosh(params.beta*updata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	    cosharraydown[k*Ny + m] = 1.0/std::cosh(params.beta*downdata[k*Ny + m].lambda/2.0)/2.0/static_cast<FPType>(params.Nx);
	  }
#pragma omp parallel for
	for(uint n = 0; n < Ny; ++n)
	{
	  FPType* tempup = new FPType[sites];
	  FPType* tempdown = new FPType[sites];
	  for(uint k = 0; k < params.Nx; ++k)
	    for(uint m = 0; m < Ny; ++m)
	    {
	      tempup[k*Ny + m] = cosharrayup[k*Ny + m] * norm(updata[k*Ny + /*n*/m].evec[/*m*/n]);
	      tempdown[k*Ny + m] = cosharraydown[k*Ny + m] * norm(downdata[k*Ny + /*n*/m].evec[/*m*/n]);
	    }
	    FPType ftup = 0.0;
	    FPType ftdown = 0.0;
	    FPType arg = - params.beta/2.0;//from that we only need a tau == 0 quantity
	    FPType cup = 0.0;
	    FPType cdown = 0.0;
	    for(uint k = 0; k < params.Nx; ++k)//being careful we employ Kahan summation
	      for(uint m = 0; m < Ny; ++m)
	      {
		FPType argup = tempup[k*Ny + m] * std::exp(updata[k*Ny + m].lambda*arg);
		FPType argdown = tempdown[k*Ny + m] * std::exp(downdata[k*Ny + m].lambda*arg);
		FPType y = argup - cup;
		FPType t = ftup + y;
		cup = (t - ftup) - y;
		ftup = t;
		y = argdown - cdown;
		t = ftdown + y;
		cdown = (t - ftdown) - y;
		ftdown = t;
	      }
	   for(uint r = 0; r < params.Nx; ++r)
	   {
	      (unmixeddata        + n*params.Nx)[r] = static_cast<float>(ftup);
              (unmixeddata + sites + n*params.Nx)[r] = static_cast<float>(ftdown);
	   }
	  delete [] tempup;
	  delete [] tempdown;
	}

	delete [] cosharrayup;
	delete [] cosharraydown;
        const unsigned int Nxw = Nw*params.Nx;
	FPType normierungunmixed = 1.0/params.beta / params.Nx;
	FPType normierungmixed = 1.0/params.beta / std::sqrt(params.Nx);
	const std::complex<FPType> expNt = std::exp(std::complex<FPType>(0.0, -M_PI/(Nt-1)));
	mixeddata = new std::complex<FPType>[2*Nxyt];
#pragma omp parallel for
        for (uint n = 0; n < Ny; ++n)
        {
	    std::complex<FPType>* Qup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* Qdown = new std::complex<FPType>[Nxw];
	    std::complex<FPType>* funup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* fundown = new std::complex<FPType>[Nxw];
            memset(funup, 0, Nxw*sizeof(std::complex<FPType>));
            memset(fundown, 0, Nxw*sizeof(std::complex<FPType>));
//            double start = omp_get_wtime();
//            std::cout<<"n = "<<n<<std::endl;
            for (uint k = 0; k < params.Nx; ++k)
            {
                for (uint m = 0; m < Ny; ++m)
                {
                    std::complex<FPType> facup   = conj(updata  [k*Ny + m].u) * updata[k*Ny + m].evec[n];
                    std::complex<FPType> facdown = conj(downdata[k*Ny + m].u) * downdata[k*Ny + m].evec[n];
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {//the layout of the frequencies is now (w_n, -w_n) , that is every negative frequency is stored next to its positive counterpart. 
                     //Hopefully this gives a better data locality
                        funup[2*omega_idx*params.Nx + k]         += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        funup[(2*omega_idx + 1)*params.Nx + k]   += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, omega[omega_idx].omega);
                        fundown[2*omega_idx*params.Nx + k]       += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        fundown[(2*omega_idx + 1)*params.Nx + k] += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, omega[omega_idx].omega);
                    }
                }
            }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint w = 0; w < Nw; ++w)
            {
                fourier1(reinterpret_cast<FPType*>(funup + w*params.Nx), params.Nx, 1);
                fourier1(reinterpret_cast<FPType*>(fundown + w*params.Nx), params.Nx, 1);
		//funup as well as fundown now contain Q(r, i omega) for a particular value of the orbital n
                for (uint r = 0; r < params.Nx; ++r)
                {
		    funup[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
		    fundown[w*params.Nx + r] *= pref; // == Q. pref == V/sqrt(L)
                    (Qup + r*Nw)[w] = funup[w*params.Nx + r];//norm(funup[w*params.Nx + r]);
                    (Qdown + r*Nw)[w] = fundown[w*params.Nx + r];//norm(fundown[w*params.Nx + r]);
//		    funup[w*params.Nx + r] = funup[w*params.Nx + r];
//		    fundown[w*params.Nx + r] = fundown[w*params.Nx + r];
                }
            }
            for(uint r = 0; r < params.Nx; ++r)
            for(uint w = 0; w < Nw/2; ++w)
	    {
	      std::complex<FPType> temp = conj((Qup + r*Nw)[2*w]);
	      (Qup + r*Nw)[2*w] *= conj((Qup + r*Nw)[2*w+1]);
	      (Qup + r*Nw)[2*w+1] *= temp;
	      temp = conj((Qdown + r*Nw)[2*w]);
	      (Qdown + r*Nw)[2*w] *= conj((Qdown + r*Nw)[2*w+1]);
	      (Qdown + r*Nw)[2*w+1] *= temp;
	    }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint r = 0; r < params.Nx; ++r)
	    {
	      std::complex<FPType> expt = 1;
	      {//an empty block for the unmixeddata
                    std::complex<FPType> tempup = 0;
                    std::complex<FPType> tempdown = 0;
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {
                        std::complex<FPType> gupp = gup[2*omega_idx];
			std::complex<FPType> gupm = gup[2*omega_idx + 1];
                        std::complex<FPType> gdownp = gdown[2*omega_idx];
			std::complex<FPType> gdownm = gdown[2*omega_idx + 1];
                        tempup += (Qup + r * Nw)[2*omega_idx] * gupp + (Qup + r * Nw)[2*omega_idx + 1] * gupm;
                        tempdown += (Qdown + r * Nw)[2*omega_idx] * gdownp + (Qdown + r * Nw)[2*omega_idx + 1] * gdownm;
                    }
                    (unmixeddata        + n*params.Nx)[r] += tempup*normierungunmixed;
                    (unmixeddata + sites + n*params.Nx)[r] += tempdown*normierungunmixed;
	      }
	      for(uint t = 0; t < Nt; ++t)// here is the final Matsubara transform
	      {
                    std::complex<FPType> tempupmixed = 0;
                    std::complex<FPType> tempdownmixed = 0;
		    std::complex<FPType> expiom = expt;
		    std::complex<FPType> expfac = expiom*expiom;
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {
                        std::complex<FPType> gupp = gup[2*omega_idx];
			std::complex<FPType> gupm = gup[2*omega_idx + 1];
                        std::complex<FPType> gdownp = gdown[2*omega_idx];
			std::complex<FPType> gdownm = gdown[2*omega_idx + 1];
			std::complex<FPType> cexpiom = conj(expiom);
			tempupmixed += cexpiom * funup[(2*omega_idx)*params.Nx + r] * gupp + expiom * funup[(2*omega_idx + 1)*params.Nx + r] * gupm;
			tempdownmixed += cexpiom * fundown[(2*omega_idx)*params.Nx + r] * gdownp + expiom * fundown[(2*omega_idx + 1)*params.Nx + r] * gdownm;
			expiom *= expfac;
                    }
                    (mixeddata        + n*params.Nx*Nt + r*Nt)[t] = tempupmixed*normierungmixed;
                    (mixeddata + Nxyt + n*params.Nx*Nt + r*Nt)[t] = tempdownmixed*normierungmixed;
		    expt *= expNt;
	      }
//                if (r > 3 )exit(-1);
//		test<<"&"<<std::endl;
	    }
//        std::cout<<"Initialization took "<<omp_get_wtime() - start<<" seconds"<<std::endl;
        delete [] funup;
        delete [] fundown;
	delete [] Qup;
	delete [] Qdown;
        }
#ifdef _OPENMP
        std::cout<<"Initialization took "<<omp_get_wtime() - start2<<" seconds"<<std::endl;
#endif
	delete [] gup;
	delete [] gdown;
	delete [] omega;
        std::cout<<"KondoCloud_X done"<<std::endl;
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func) 
    {
    contrib_dry(UP, UP, func);
    contrib_dry(UP, DOWN, func);
    contrib_dry(DOWN, UP, func);
    contrib_dry(DOWN, DOWN, func);      
    }
    /**
    this determines the KondoCloud_X for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& functionpoints;
    const FPType beta;
    std::complex<FPType>* mixeddata;///< here we store <d^+ gamma>.
    std::complex<FPType>* unmixeddata;///< here we store <gamma^+ gamma>
    const uint Ny;
    const uint Nx;
    const uint Nt;
    const uint Nw;
    const FPType delta_beta;
    const uint32_t sites;
    const uint Nxyt;
    std::complex<FPType> accessmixed(uint r, uint n, SPINS spin1, FPType tau1 , SPINS spin2, FPType tau2) const
    {
      FPType delta_tau = tau1 - tau2;
      FPType sign = 1.0;
      std::complex<FPType>* dataptr = mixeddata;
      if(spin1 != spin2) return 0.0;//FIXME!!!!!!!!!!!!!!!!!!!! only the case for the spin symmetric case!
      if (spin1 == DOWN) dataptr += Nxyt;
      dataptr = dataptr + n*Nx*Nt + r*Nt;
//      if(std::abs(delta_tau) < std::numeric_limits<FPType>::epsilon())
      if(fpequal(tau1, tau2))
	return std::complex<FPType>(dataptr[0]);
      if(delta_tau < 0)
      {
	sign = -1.0;
	delta_tau += beta;
      }
      FPType fptau_idx0;
      FPType rem = std::modf(delta_tau/delta_beta, &fptau_idx0);//round to the smaller index and determine how far we're off.
      std::size_t tau_idx0 = lround(fptau_idx0);
      return std::complex<FPType>(lerp(rem, dataptr[tau_idx0], dataptr[tau_idx0 + 1]))*sign;
    }
    std::complex<FPType> accessunmixed(uint r, uint n, SPINS spin) const
    {
      std::complex<FPType>* dataptr = unmixeddata;
      if (spin == DOWN) dataptr += sites;
      dataptr = dataptr + n*Nx + r;
      return dataptr[0];
    }
    void contrib_dry(SPINS sigma, SPINS sigmaprime, DryRun<typename Configuration::value_type, GFRetVal>& func)
    {//yes, the only thing that we don't derive from tables only depends on sigmaprime 
        typename GreensFunction::Vertex v1;
        v1.spin = sigmaprime;
        v1.tau = 0;
        func(v1, v1);
    }
    GFRetVal dotghelper(FPType tau1, SPINS spin1, FPType tau2, SPINS spin2, const Configuration config)
    {
//    auto gdot = [](FPType tau1, SPINS spin1, FPType tau2, SPINS spin2){return GreensFunction::eval(typename GreensFunction::Vertex(tau1, spin1), typename GreensFunction::Vertex(tau2, spin2));};
    struct GDot{GFRetVal operator()
    (FPType tau1, SPINS spin1, FPType tau2, SPINS spin2){return GreensFunction::eval(typename GreensFunction::Vertex(tau1, spin1), typename GreensFunction::Vertex(tau2, spin2));}
    } gdot;
    GFRetVal retval = gdot(tau1, spin1, tau2, spin2);
    for(uint r = 0; r < config.size(); ++r)
      for(uint s = 0; s < config.size(); ++s)
      {
	retval -= gdot(config[r].tau, UP, tau1, spin1) * config.matcont.mat(2*r, 2*s) * gdot(tau2, spin2, config[s].tau, UP);
	retval -= gdot(config[r].tau, DOWN, tau1, spin1) * config.matcont.mat(2*r+1, 2*s) * gdot(tau2, spin2, config[s].tau, UP);
	retval -= gdot(config[r].tau, UP, tau1, spin1) * config.matcont.mat(2*r, 2*s+1) * gdot(tau2, spin2, config[s].tau, DOWN);
	retval -= gdot(config[r].tau, DOWN, tau1, spin1) * config.matcont.mat(2*r+1, 2*s+1) * gdot(tau2, spin2, config[s].tau, DOWN);
      }
      return retval;
    }
};

template <class Config, class GreensFunction, SPINS Spin>
void KondoCloud_X<Config, GreensFunction, Spin>::evaluate(const Configuration& config, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(sites);
    typename GreensFunction::Vertex v1(0.0, UP);
    typename GreensFunction::Vertex v2(0.0, DOWN);
    GFRetVal dotcontrib = dowick(v1, v1) - dowick(v2,v2);
    GFRetVal *const dotgf = new GFRetVal[config.size()*2*2*2];
    for(uint i = 0; i < config.size(); ++i)//here we sum up < d^+ d(t) >
    {
	dotgf[4*i + 0] = dotghelper(0, UP, config[i].tau, UP, config);
	dotgf[4*i + 1] = dotghelper(0, UP, config[i].tau, DOWN, config);
	dotgf[4*i + 2] = dotghelper(0, DOWN, config[i].tau, UP, config);
	dotgf[4*i + 3] = dotghelper(0, DOWN, config[i].tau, DOWN, config);
    }
    for(uint i = 0; i < config.size(); ++i)//here we sum up < d^+ d(t) >
    {
	dotgf[4*(i+config.size()) + 0] = dotghelper(config[i].tau, UP, 0, UP, config);
	dotgf[4*(i+config.size()) + 1] = dotghelper(config[i].tau, UP, 0, DOWN, config);
	dotgf[4*(i+config.size()) + 2] = dotghelper(config[i].tau, DOWN, 0, UP, config);
	dotgf[4*(i+config.size()) + 3] = dotghelper(config[i].tau, DOWN, 0, DOWN, config);
    }
//    std::ofstream kc("kc.txt");
#pragma omp parallel for
    for (unsigned int k = 0; k < func.size(); ++k)
    {
	uint n = k / Nx;
	uint r = k % Nx;
// 	GFRetVal unmixedup = accessunmixed(r, n, UP);
// 	GFRetVal unmixeddown = accessunmixed(r, n, DOWN);
	GFRetVal mixedupup =  accessmixed(r, n, UP, 0.0, UP, 0.0);
// 	GFRetVal mixedupdown = accessmixed(r, n, UP, 0.0, DOWN, 0.0);
// 	GFRetVal mixeddownup = accessmixed(r, n, DOWN, 0.0, UP, 0.0);
	GFRetVal mixeddowndown = accessmixed(r, n, DOWN, 0.0, DOWN, 0.0);

	GFRetVal mixedupupconj =  conj(accessmixed(r, n, UP, 0.0, UP, 0.0));
// 	GFRetVal mixedupdownconj = conj(accessmixed(r, n, UP, 0.0, DOWN, 0.0));
// 	GFRetVal mixeddownupconj = conj(accessmixed(r, n, DOWN, 0.0, UP, 0.0));
	GFRetVal mixeddowndownconj = conj(accessmixed(r, n, DOWN, 0.0, DOWN, 0.0));
	for(uint q = 0; q < config.size(); ++q)
	  for(uint s = 0; s < config.size(); ++s)
	  {
// 	    unmixedup -= (
// 	      config.matcont.mat(2*q, 2*s)  *    accessmixed(r, n, UP, config[q].tau, UP, 0.0)  *conj(accessmixed(r, n, UP, config[s].tau, UP, 0.0))
// 	     +config.matcont.mat(2*q, 2*s+1)*    accessmixed(r, n, UP, config[q].tau, UP, 0.0)  *conj(accessmixed(r, n, UP, config[s].tau, DOWN, 0.0))
// 	     +config.matcont.mat(2*q+1, 2*s)*    accessmixed(r, n, DOWN, config[q].tau, UP, 0.0)*conj(accessmixed(r, n, UP, config[s].tau, UP, 0.0))
// 	     +config.matcont.mat(2*q+1, 2*s+1)  *accessmixed(r, n, DOWN, config[q].tau, UP, 0.0)*conj(accessmixed(r, n, UP, config[s].tau, DOWN, 0.0))
// 	    );
// 	    unmixeddown -= (
// 	      config.matcont.mat(2*q, 2*s)  *    accessmixed(r, n, UP, config[q].tau, DOWN, 0.0)  *conj(accessmixed(r, n, DOWN, config[s].tau, UP, 0.0))
// 	     +config.matcont.mat(2*q, 2*s+1)*    accessmixed(r, n, UP, config[q].tau, DOWN, 0.0)  *conj(accessmixed(r, n, DOWN, config[s].tau, DOWN, 0.0))
// 	     +config.matcont.mat(2*q+1, 2*s)*    accessmixed(r, n, DOWN, config[q].tau, DOWN, 0.0)*conj(accessmixed(r, n, DOWN, config[s].tau, UP, 0.0))
// 	     +config.matcont.mat(2*q+1, 2*s+1)  *accessmixed(r, n, DOWN, config[q].tau, DOWN, 0.0)*conj(accessmixed(r, n, DOWN, config[s].tau, DOWN, 0.0))
// 	    );
	    
	    typename GreensFunction::Vertex v1u(0.0, UP);
	    typename GreensFunction::Vertex v2su(config[s].tau, UP);
	    typename GreensFunction::Vertex v1d(0.0, DOWN);
	    typename GreensFunction::Vertex v2sd(config[s].tau, DOWN);
	    mixedupup -= (
	      config.matcont.mat(2*q, 2*s)  *    accessmixed(r, n, UP, config[q].tau, UP, 0.0)  * dotgf[4*s+0]
	     +config.matcont.mat(2*q, 2*s+1)*    accessmixed(r, n, UP, config[q].tau, UP, 0.0)  * dotgf[4*s+1]
	     +config.matcont.mat(2*q+1, 2*s)*    accessmixed(r, n, DOWN, config[q].tau, UP, 0.0)* dotgf[4*s+0]
	     +config.matcont.mat(2*q+1, 2*s+1)  *accessmixed(r, n, DOWN, config[q].tau, UP, 0.0)* dotgf[4*s+1]
	    );

// 	    mixedupdown -= (
// 	      config.matcont.mat(2*q, 2*s)  *    accessmixed(r, n, UP, config[q].tau, DOWN, 0.0)  * dotgf[4*s+0]
// 	     +config.matcont.mat(2*q, 2*s+1)*    accessmixed(r, n, UP, config[q].tau, DOWN, 0.0)  * dotgf[4*s+1]
// 	     +config.matcont.mat(2*q+1, 2*s)*    accessmixed(r, n, DOWN, config[q].tau, DOWN, 0.0)* dotgf[4*s+0]
// 	     +config.matcont.mat(2*q+1, 2*s+1)  *accessmixed(r, n, DOWN, config[q].tau, DOWN, 0.0)* dotgf[4*s+1]
// 	    );
// 	    
// 	     mixeddownup -= (
// 	      config.matcont.mat(2*q, 2*s)  *    accessmixed(r, n, UP, config[q].tau, UP, 0.0)  * dotgf[4*s+2]
// 	     +config.matcont.mat(2*q, 2*s+1)*    accessmixed(r, n, UP, config[q].tau, UP, 0.0)  * dotgf[4*s+3]
// 	     +config.matcont.mat(2*q+1, 2*s)*    accessmixed(r, n, DOWN, config[q].tau, UP, 0.0)* dotgf[4*s+2]
// 	     +config.matcont.mat(2*q+1, 2*s+1)  *accessmixed(r, n, DOWN, config[q].tau, UP, 0.0)* dotgf[4*s+3]
// 	    );

   	    mixeddowndown -= (
	      config.matcont.mat(2*q, 2*s)  *    accessmixed(r, n, UP, config[q].tau, DOWN, 0.0)  * dotgf[4*s+2]
	     +config.matcont.mat(2*q, 2*s+1)*    accessmixed(r, n, UP, config[q].tau, DOWN, 0.0)  * dotgf[4*s+3]
	     +config.matcont.mat(2*q+1, 2*s)*    accessmixed(r, n, DOWN, config[q].tau, DOWN, 0.0)* dotgf[4*s+2]
	     +config.matcont.mat(2*q+1, 2*s+1)  *accessmixed(r, n, DOWN, config[q].tau, DOWN, 0.0)* dotgf[4*s+3]
	    );
	    
	    typename GreensFunction::Vertex v2qu(config[q].tau, UP);
	    typename GreensFunction::Vertex v2qd(config[q].tau, DOWN);
	    mixedupupconj -= (
	      config.matcont.mat(2*q, 2*s)  *    conj(accessmixed(r, n, UP, 0.0, UP, config[s].tau))  * dotgf[4*(q+config.size()) + 0]
	     +config.matcont.mat(2*q, 2*s+1)*    conj(accessmixed(r, n, UP, 0.0, DOWN, config[s].tau))  * dotgf[4*(q+config.size()) + 0]
	     +config.matcont.mat(2*q+1, 2*s)*    conj(accessmixed(r, n, UP, 0.0, UP, config[s].tau))* dotgf[4*(q+config.size()) + 2]
	     +config.matcont.mat(2*q+1, 2*s+1)  *conj(accessmixed(r, n, UP, 0.0, DOWN, config[s].tau))* dotgf[4*(q+config.size()) + 2]
	    );

// 	    mixedupdownconj -= (
// 	      config.matcont.mat(2*q, 2*s)  *    conj(accessmixed(r, n, UP, 0.0, UP, config[s].tau))  * dotgf[4*(q+config.size()) + 1]
// 	     +config.matcont.mat(2*q, 2*s+1)*    conj(accessmixed(r, n, UP, 0.0, DOWN, config[s].tau))  * dotgf[4*(q+config.size()) + 1]
// 	     +config.matcont.mat(2*q+1, 2*s)*    conj(accessmixed(r, n, UP, 0.0, UP, config[s].tau))* dotgf[4*(q+config.size()) + 3]
// 	     +config.matcont.mat(2*q+1, 2*s+1)  *conj(accessmixed(r, n, UP, 0.0, DOWN, config[s].tau))* dotgf[4*(q+config.size()) + 3]
// 	    );
// 	    
// 	     mixeddownupconj -= (
// 	      config.matcont.mat(2*q, 2*s)  *    conj(accessmixed(r, n, DOWN, 0.0, UP, config[s].tau))  * dotgf[4*(q+config.size()) + 0]
// 	     +config.matcont.mat(2*q, 2*s+1)*    conj(accessmixed(r, n, DOWN, 0.0, DOWN, config[s].tau))  * dotgf[4*(q+config.size()) + 0]
// 	     +config.matcont.mat(2*q+1, 2*s)*    conj(accessmixed(r, n, DOWN, 0.0, UP, config[s].tau))* dotgf[4*(q+config.size()) + 2]
// 	     +config.matcont.mat(2*q+1, 2*s+1)  *conj(accessmixed(r, n, DOWN, 0.0, DOWN, config[s].tau))* dotgf[4*(q+config.size()) + 2]
// 	    );

   	    mixeddowndownconj -= (
	      config.matcont.mat(2*q, 2*s)  *    conj(accessmixed(r, n, DOWN, 0.0, UP, config[s].tau))  * dotgf[4*(q+config.size()) + 1]
	     +config.matcont.mat(2*q, 2*s+1)*    conj(accessmixed(r, n, DOWN, 0.0, DOWN, config[s].tau))  * dotgf[4*(q+config.size()) + 1]
	     +config.matcont.mat(2*q+1, 2*s)*    conj(accessmixed(r, n, DOWN, 0.0, UP, config[s].tau))* dotgf[4*(q+config.size()) + 3]
	     +config.matcont.mat(2*q+1, 2*s+1)  *conj(accessmixed(r, n, DOWN, 0.0, DOWN, config[s].tau))* dotgf[4*(q+config.size()) + 3]
	    );
	  }
            //add to measurement
            //note that this formula is the relevant one if the dot greens function and the mixed <dc> greensfunctions are spin diagonal
            func[k] = ( - mixedupup*mixeddowndownconj - mixeddowndown*mixedupupconj)* config.phase/4.0;
//kc<<dotcontrib <<" "<< (unmixedup - unmixeddown)<<" "<<mixedupup<<" "<<mixedupupconj<<" "<< mixedupdown<<" "<<mixedupdownconj <<" "<< mixeddownup<<" "<<mixeddownupconj <<" "<< mixeddowndown<<" "<<mixeddowndownconj<<std::endl;
    }
//    exit(-1);
    this->add_bin(func);
    delete [] dotgf;
    return;
}

/**
A class for measuring the hybridization:
<d^\dagger_{-\sigma} a_{0,-sigma}>
*/
template <class Config, class GreensFunction, SPINS Spin>
class Hybridization : public Network_Cache<Config, typename Config::SignType>
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef GreensFunction GF;
    typedef std::valarray<GFRetVal> Function;
    typedef GFRetVal ObservableType;///<the Hybridization is only a complex number
    /**
    The Constructor for the Hybrization
    */
    Hybridization(typename Config::Comm& net, const Parameters& params)/* throw()*/ : Network_Cache<Config, ObservableType>(net, "Hybridization"), beta(params.beta), Ny(params.Nb*2), Nx(params.Nx), Nt(400 + 1), Nw(4*Nt),
            delta_beta(params.beta/static_cast<FPType>(Nt-1))
    {
#ifdef _OPENMP
        double start2 = omp_get_wtime();
#endif
        std::cout<<"Creating Data for Hybridization"<<std::endl;
        const GOmegaData<FPType> *const updata = GF::gomegaup->data;
        const GOmegaData<FPType> *const downdata = GF::gomegadown->data;
        Omegadata<FPType>* omega = new Omegadata<FPType>[Nw/2];
	std::complex<FPType>* gup = new std::complex<FPType>[Nw];
	std::complex<FPType>* gdown = new std::complex<FPType>[Nw];
        for (int t = 0; t < static_cast<int>(Nw/2); ++t)
	{
	    gup[2*t] = conj((*GF::gomegaup)(t));
	    gup[2*t + 1] = conj((*GF::gomegaup)( -t - 1));
	    gdown[2*t] = conj((*GF::gomegadown)(t));
	    gdown[2*t + 1] = conj((*GF::gomegadown)( -t - 1));
            omega[t].omega = M_PI/params.beta*(2*t + 1);
	    omega[t].omegasq = omega[t].omega*omega[t].omega;
	    omega[t].invomega = 1.0/omega[t].omega;
	}
        FPType pref = GF::gethybridization() / std::sqrt(params.Nx);
        const unsigned int Nxw = Nw*params.Nx;
	FPType normierungmixed = 1.0/params.beta / std::sqrt(params.Nx);
	const std::complex<FPType> expNt = std::exp(std::complex<FPType>(0.0, -M_PI/(Nt-1)));
	mixeddata = new std::complex<FPType>[2*Nt];
	    std::complex<FPType>* funup = new std::complex<FPType>[Nxw];
            std::complex<FPType>* fundown = new std::complex<FPType>[Nxw];
            memset(funup, 0, Nxw*sizeof(std::complex<FPType>));
            memset(fundown, 0, Nxw*sizeof(std::complex<FPType>));
//            double start = omp_get_wtime();
            for (uint k = 0; k < params.Nx; ++k)
            {
                for (uint m = 0; m < Ny; ++m)
                {
                    std::complex<FPType> facup   = conj(updata  [k*Ny + m].u) * updata[k*Ny + m].evec[0];
                    std::complex<FPType> facdown = conj(downdata[k*Ny + m].u) * downdata[k*Ny + m].evec[0];
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {//the layout of the frequencies is now (w_n, -w_n) , that is every negative frequency is stored next to its positive counterpart. 
                     //Hopefully this gives a better data locality
                        funup[2*omega_idx*params.Nx + k]         += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        funup[(2*omega_idx + 1)*params.Nx + k]   += facup/std::complex<FPType>(-updata[k*Ny + m].lambda, omega[omega_idx].omega);
                        fundown[2*omega_idx*params.Nx + k]       += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, -omega[omega_idx].omega);
                        fundown[(2*omega_idx + 1)*params.Nx + k] += facdown/std::complex<FPType>(-downdata[k*Ny + m].lambda, omega[omega_idx].omega);
                    }
                }
            }
//            std::cout<<"time now: "<<omp_get_wtime() - start<<std::endl;
            for (uint w = 0; w < Nw; ++w)
            {
                fourier1(reinterpret_cast<FPType*>(funup + w*params.Nx), params.Nx, 1);
                fourier1(reinterpret_cast<FPType*>(fundown + w*params.Nx), params.Nx, 1);
		//funup as well as fundown now contain Q(r, i omega) for for n = 0
		    funup[w*params.Nx] *= pref; // == Q. pref == V/sqrt(L)
		    fundown[w*params.Nx] *= pref; // == Q. pref == V/sqrt(L)
            }
	      std::complex<FPType> expt = 1;
	      for(uint t = 0; t < Nt; ++t)// here is the final Matsubara transform
	      {
                    std::complex<FPType> tempupmixed = 0;
                    std::complex<FPType> tempdownmixed = 0;
		    std::complex<FPType> expiom = expt;
		    std::complex<FPType> expfac = expiom*expiom;
                    for (int omega_idx = 0; omega_idx < static_cast<int>(Nw/2); ++omega_idx)
                    {
                        std::complex<FPType> gupp = gup[2*omega_idx];
			std::complex<FPType> gupm = gup[2*omega_idx + 1];
                        std::complex<FPType> gdownp = gdown[2*omega_idx];
			std::complex<FPType> gdownm = gdown[2*omega_idx + 1];
			std::complex<FPType> cexpiom = conj(expiom);
			tempupmixed += cexpiom * funup[(2*omega_idx)*params.Nx] * gupp + expiom * funup[(2*omega_idx + 1)*params.Nx] * gupm;
			tempdownmixed += cexpiom * fundown[(2*omega_idx)*params.Nx] * gdownp + expiom * fundown[(2*omega_idx + 1)*params.Nx] * gdownm;
			expiom *= expfac;
                    }
                    mixeddata[t] = tempupmixed*normierungmixed;
                    (mixeddata + Nt)[t] = tempdownmixed*normierungmixed;
		    expt *= expNt;
	      }
//        std::cout<<"Initialization took "<<omp_get_wtime() - start<<" seconds"<<std::endl;
        delete [] funup;
        delete [] fundown;
#ifdef _OPENMP
        std::cout<<"Initialization took "<<omp_get_wtime() - start2<<" seconds"<<std::endl;
#endif
	std::ofstream file("goff.txt");
	for(uint k = 0; k < Nt; ++k)
	  file<<mixeddata[k]<<std::endl;
	delete [] gup;
	delete [] gdown;
	delete [] omega;
        std::cout<<"Hybridization done"<<std::endl;
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func) 
    {
    }
    /**
    this determines the Hybridization for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const FPType beta;
    std::complex<FPType>* mixeddata;///< here we store <d^+ gamma>.
    const uint Ny;
    const uint Nx;
    const uint Nt;
    const uint Nw;
    const FPType delta_beta;
    std::complex<FPType> accessmixed(SPINS spin1, FPType tau1 , SPINS spin2, FPType tau2) const
    {
      FPType delta_tau = tau1 - tau2;
      FPType sign = 1.0;
      std::complex<FPType>* dataptr = mixeddata;
      if(spin1 != spin2) return 0.0;//FIXME!!!!!!!!!!!!!!!!!!!! only the case for the spin symmetric case!
      if (spin1 == DOWN) dataptr += Nt;
      if(std::abs(delta_tau) < std::numeric_limits<FPType>::epsilon())
	return std::complex<FPType>(dataptr[0]);
      if(delta_tau < 0)
      {
	sign = -1.0;
	delta_tau += beta;
      }
      FPType fptau_idx0;
      FPType rem = std::modf(delta_tau/delta_beta, &fptau_idx0);//round to the smaller index and determine how far we're off.
      std::size_t tau_idx0 = lround(fptau_idx0);
      return std::complex<FPType>(lerp(rem, dataptr[tau_idx0], dataptr[tau_idx0 + 1]))*sign;
    }
};

template <class Config, class GreensFunction, SPINS Spin>
void Hybridization<Config, GreensFunction, Spin>::evaluate(const Configuration& config, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType retval = accessmixed(!Spin, 0, !Spin, 0);
	struct GDot{GFRetVal operator()
    (FPType tau2, SPINS spin2){return GreensFunction::eval(typename GreensFunction::Vertex(0, !Spin), typename GreensFunction::Vertex(tau2, spin2));}
    } gdot;
	for(uint q = 0; q < config.size(); ++q)
	  for(uint s = 0; s < config.size(); ++s)
	  {
	    retval -= (
	      config.matcont.mat(2*q, 2*s)  *  accessmixed(UP, config[q].tau, !Spin, 0.0)  * gdot(config[s].tau, UP)
	     +config.matcont.mat(2*q, 2*s+1)*  accessmixed(UP, config[q].tau, !Spin, 0.0)  * gdot(config[s].tau, DOWN)
	     +config.matcont.mat(2*q+1, 2*s)*  accessmixed(DOWN, config[q].tau, !Spin, 0.0)* gdot(config[s].tau, UP)
	     +config.matcont.mat(2*q+1, 2*s+1)*accessmixed(DOWN, config[q].tau, !Spin, 0.0)* gdot(config[s].tau, DOWN)
	    );
	  }
    //add to measurement
    this->add_bin(retval*config.phase);
    return;
}

/**
 * Depending on s this measures the charge charge( s=1) correlation function or the spin-spin (s=-1) correlation function
 * Note also that by its bare definition it is a real quantity.
 * */
template <class Config, int s>
class SpinChargeParent : public Network_Cache<Config, std::valarray<std::valarray<typename Config::Configuration::FPType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<FPType> Function;
    typedef std::valarray<Function> ObservableType;///< spin and charge correlations are potentially complex in kspace and are spatially resolved time-dependent observable
    /**
    The Constructor
    */
    SpinChargeParent(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config,ObservableType>(n, (s==1)?"CharcheChargeCorrelation":"SpinSpinCorrelation"),
    len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the observable for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config, int cs>
void SpinChargeParent<Config, cs>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    const typename Configuration::value_type v2(0, 0, UP);
    const typename Configuration::value_type v4(0, 0, DOWN);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            FPType s = j * delta_s;
	    if(j == 0)
	      s = 0.000001;//seems to be necessary here...
            for(uint r = 0; r < len; ++r)
            {
		const typename Configuration::value_type v1(r, s, UP);
		const typename Configuration::value_type v3(r, s, DOWN);
		genericTwoParticleGreensfunction_dry<Configuration>(func, v1, v1, v2, v2);
		genericTwoParticleGreensfunction_dry<Configuration>(func, v3, v3, v4, v4); 
		genericTwoParticleGreensfunction_dry<Configuration>(func, v1, v1, v4, v4); 
		genericTwoParticleGreensfunction_dry<Configuration>(func, v3, v3, v2, v2);
            }
        }
    return;
}

template <typename T>
struct Help
{
  typedef T RetType;
static inline RetType toreal(T a)
{
  return a;
}
};

template <typename FPType>
struct Help<std::complex<FPType> >
{
  typedef FPType RetType;
static inline RetType toreal(std::complex<FPType> a)
{
  return a.real();
}
};

template <class Config, int cs>
void SpinChargeParent<Config, cs>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    const FPType fac = static_cast<FPType>(TWOPI / len);
    const FPType csfac = (cs == 1? 1.0: 0.25);//switch prefactors
    const typename Config::SignType phase_factor = configuration.phase * csfac;
    const typename Configuration::value_type v2(0, 0, UP);
    const typename Configuration::value_type v4(0, 0, DOWN);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            FPType s = j * delta_s;
	    if(j == 0)
	      s = 0.000001;//seems to be necessary here, to get the correct function value at tau=0...
	    GFRetVal sum = 0;
            for(uint r = 0; r < len; ++r)
            {
                FPType pref = std::cos(fac * (k*r));
		const typename Configuration::value_type v1(r, s, UP);
		const typename Configuration::value_type v3(r, s, DOWN);
		auto sum2 = genericTwoParticleGreensfunction<Configuration>(dowick, v1, v1, v2, v2);//to deduce type
		sum2 += genericTwoParticleGreensfunction<Configuration>(dowick, v3, v3, v4, v4) 
		+ FPType(cs)*(genericTwoParticleGreensfunction<Configuration>(dowick, v1, v1, v4, v4) 
		   + genericTwoParticleGreensfunction<Configuration>(dowick, v3, v3, v2, v2));
                sum += pref * sum2;
            }
            func[k][j] = Help<GFRetVal>::toreal(phase_factor * sum);//normalization not necessary since the 1/N factor is already in the definition of G_0
        }
    }
    this->add_bin(func);
    return;
}

/**
 * This measures <S^+(k, tau) S^-(0,0)> which is due to translation symmetry the same as <S^+(k,tau), S^-(k,0)>
 * Note also that by its bare definition it is a real quantity.
 * */
template <class Config>
class SplusSminus : public Network_Cache<Config, std::valarray<std::valarray<typename Config::Configuration::FPType> > >
{
public:
    typedef typename Config::Configuration Configuration;
    typedef typename Configuration::FPType FPType;
    typedef typename Config::SignType GFRetVal;
    typedef std::valarray<FPType> Function;///< SplusSminus is a real quantity
    typedef std::valarray<Function> ObservableType;///< spin and charge correlations are potentially complex in kspace and are spatially resolved time-dependent observable
    /**
    The Constructor
    */
    SplusSminus(typename Config::Comm& n, const Parameters& params) throw() : Network_Cache<Config,ObservableType>(n, "SplusSminus"),
    len(params.N), functionpoints(params.functionpoints), delta_s(params.delta_s)
    {
    }
    void dryrun(DryRun<typename Configuration::value_type, GFRetVal>&);
    /**
    this determines the observable for a given order
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&, const DoWick<typename Configuration::value_type, GFRetVal>&);
private:
    const uint32_t& len;
    const uint32_t& functionpoints;
    const double delta_s;
};

template <class Config>
void SplusSminus<Config>::dryrun(DryRun<typename Configuration::value_type, GFRetVal>& func)
{
    const typename Configuration::value_type v2(0, 0, UP);
    const typename Configuration::value_type v4(0, 0, DOWN);
//    func(v4, v2);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
            FPType s = j * delta_s;
	    	    if(j == 0)
	      s = 0.000001;
            for(uint r = 0; r < len; ++r)
            {
		const typename Configuration::value_type v1(r, s, UP);
		const typename Configuration::value_type v3(r, s, DOWN);
//		func.template onSector<UP>(r, s, 0,0);
//		func.template onSector<DOWN>(0, 0, r, s);
/*func(v1, v3);
func(v1, v2);
func(v4, v3);*/
genericTwoParticleGreensfunction_dry<Configuration>(func, v1, v3, v4, v2);
            }
        }
    return;
}

template <class Config>
void SplusSminus<Config>::evaluate(const Configuration& configuration, const DoWick<typename Configuration::value_type, GFRetVal>& dowick)
{
    ObservableType func(len);
    FPType invlen = 1.0/len;
    const FPType fac = TWOPI * invlen;
    const typename Configuration::value_type v2(0, 0, UP);
    const typename Configuration::value_type v4(0, 0, DOWN);
//    auto constres = dowick(v4, v2);
    for (unsigned int k = 0; k < len; ++k)
    {
        func[k].resize(functionpoints);
        for (unsigned int j = 0; j < functionpoints; ++j)
        {
             FPType s = j * delta_s;
	    if(j == 0)
	      s = 0.000001;//seems to be necessary here...
	    GFRetVal sum = 0;
            for(uint r = 0; r < len; ++r)
            {
                FPType pref = std::cos(fac * (k*r));
		const typename Configuration::value_type v1(r, s, UP);
		const typename Configuration::value_type v3(r, s, DOWN);
		auto sum2 = 
//		dowick(v1, v3) * constres - dowick(v1, v2) * dowick(v4, v3);
//dowick.template onSector<UP>(r, s, 0,0) * dowick.template onSector<DOWN>(0, 0, r, s);
genericTwoParticleGreensfunction<Configuration>(dowick, v1, v3, v4, v2);
                sum += pref * sum2;
            }
            func[k][j] = Help<GFRetVal>::toreal(configuration.phase * sum);//normalization not necessary since the 1/N factor is already in the definition of G_0
        }
    }
    this->add_bin(func);
    return;
}
#endif
