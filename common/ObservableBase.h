/***************************************************************************
 *   Copyright (C) 2009 - 2011 by Florian Goth   *
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
#ifndef OBSERVABLEBASE_H
#define OBSERVABLEBASE_H
#include <string>
#include <complex>
#include <vector>
#include <valarray>
#include "CacheGlue.h"

template <class Obs>
struct Mean
{
    static Obs mean(const std::vector<Obs>& cont) MTL_PURE_FUNCTION
    {
        Obs ret(cont[0]);
        for (unsigned int k = 1; k < cont.size(); ++k)
            ret += cont[k];
        return ret / Obs(cont.size());
    }
};

template <class Obs>
struct Mean<std::valarray<Obs> >
{
    static std::valarray<Obs> mean(const std::vector<std::valarray<Obs> >& cont) MTL_PURE_FUNCTION
    {
        std::valarray<Obs> ret(cont[0]);
        for (unsigned int k = 1; k < cont.size(); ++k)
            ret += cont[k];
        return ret / Obs(cont.size());
    }
};

template <class Obs>
struct Mean<std::valarray<std::valarray<Obs> > >
{
    static std::valarray<std::valarray<Obs> > mean(const std::vector<std::valarray<std::valarray<Obs> > >& cont) MTL_PURE_FUNCTION
    {
        std::valarray<std::valarray<Obs> > ret(cont[0]);
        for (unsigned int k = 1; k < cont.size(); ++k)
            ret += cont[k];
        for (unsigned int i = 0; i < ret.size(); ++i)
            for (unsigned int j = 0; j < ret[i].size(); ++j)
                ret[i][j] /= Obs(cont.size());
        return ret;
    }
};

/**
Abstract base class for all Observables
*/
template <class C, class Sign>
class ObservableBase
{
public:
    typedef C Configuration;///< a typedef for the type of the Configuration
    /**
    virtual Destructor. For correctly tidying up all derived objects
    */
    ObservableBase(std::string b) : name(b) {}
    virtual ~ObservableBase() {}
    virtual void dryrun(DryRun<typename C::value_type, Sign>&) = 0;
    virtual void evaluate(const Configuration&, const DoWick<typename C::value_type, Sign>&) = 0;
    virtual void senddata() = 0;
    const std::string name;
    typedef Sign GFRetVal;
protected:
private:
};

template <class Obs_Config, class Obs_Type>
class Network_Cache : public ObservableBase<typename Obs_Config::Configuration, typename Obs_Config::SignType>
{
public:
    typedef typename Obs_Config::Comm Net;
    inline void senddata();
    inline void add_bin( const Obs_Type& val);
    inline Network_Cache(Net& n, std::string name) : ObservableBase<typename Obs_Config::Configuration, typename Obs_Config::SignType>(name), net(n) {}
    virtual ~Network_Cache() {}
protected:
    std::vector<Obs_Type> bincache;
    Net& net;
private:
};

template <class Obs_Config, class Obs_Type>
void Network_Cache<Obs_Config, Obs_Type>::add_bin(const Obs_Type& val)
{
    bincache.push_back(val);
}

template <class Obs_Config, class Obs_Type>
void Network_Cache<Obs_Config, Obs_Type>::senddata()
{
    Obs_Type temp (Mean<Obs_Type>::mean(bincache));
    net.template sendtoserver<Obs_Type>(temp);
    bincache.clear();
    return;
}
#endif
