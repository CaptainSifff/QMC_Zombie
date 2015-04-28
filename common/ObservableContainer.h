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
#ifndef OBSERVABLECONTAINER_H
#define OBSERVABLECONTAINER_H
#include "ObservableBase.h"
#include "CacheGlue.h"
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include <valarray>

template <class, class>
class ObservableBase;//forward Declaration of the base class

template<class Configuration, class GreensFunction, typename FPType>
class ObservableContainer
{
public:
    typedef typename Configuration::value_type Vertex;
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType GFRetVal;
    /**
    Add an Observable to the Container with the observables
    @param arg the Observable that gets added
    */
    void add_Observable(ObservableBase<Configuration, GFRetVal>* arg);
    void add_Cache_Value(Vertex& v1, Vertex& v2);
    /**
    This gives access to the observables
    @param k the index of the observable
    @return a pointer to the observable to allow for polymorphic behaviour
    */
    ObservableBase<Configuration, GFRetVal>*& operator[](unsigned int k) const;
    unsigned int num_Observables() const throw() MTL_PURE_FUNCTION;
    typename std::vector<ObservableBase<Configuration, GFRetVal>*>::size_type num_Cache_Entries() const throw() MTL_PURE_FUNCTION;
    /**
    This function performs the dry run with each observable and with OpenMP resizes the items storage
    */
    void dryrun();
    /**
    This Function does the actual measurement process. It evaluates all needed Greensfunction-entries and executes then the measurement procedures of the stored observables
    @param conf the Configuration for which to evaluate the Greensfunction
    */
    inline void measure(const Configuration& conf);
    void dump_cache() const;
    inline void executefunction(void (ObservableBase<Configuration, GFRetVal>::*f)());
    ~ObservableContainer();
    template <class CFG>
    ObservableContainer(const CFG& cfg);
private:
    Cache<Vertex, GFRetVal> cache;///< the cache with all the requested greensfunction values
    std::vector<ObservableBase<Configuration, GFRetVal>*> container;///< the container that contains all the observables we have registered
};

template <class Configuration, class GreensFunction, typename FPType>
template <class CFG>
ObservableContainer<Configuration, GreensFunction, FPType>::ObservableContainer(const CFG& cfg) : cache(cfg.N)
{
}

template <class Configuration, class GreensFunction, typename FPType>
ObservableContainer<Configuration, GreensFunction, FPType>::~ObservableContainer()
{
    for (std::size_t k = 0; k < container.size(); ++k)
        delete container[k];
    return;
}

template <class Configuration, class GreensFunction, typename FPType>
void ObservableContainer<Configuration, GreensFunction, FPType>::executefunction(void (ObservableBase<Configuration, GFRetVal>::*f)())
{
    for_each(container.begin(), container.end(), mem_fun(f));
//(container[0]->*f)(); how to use member to function pointers
}

template <class Configuration, class GreensFunction, typename FPType>
typename std::vector<ObservableBase<Configuration, typename GreensFunction::FreeGreensFunctionReturnValueType>*>::size_type ObservableContainer<Configuration, GreensFunction, FPType>::num_Cache_Entries() const throw()
{
    return cache.size();
}

template <class Configuration, class GreensFunction, typename FPType>
unsigned int ObservableContainer<Configuration, GreensFunction, FPType>::num_Observables() const throw()
{
    return container.size();
}

template <class Configuration, class GreensFunction, typename FPType>
void ObservableContainer<Configuration, GreensFunction, FPType>::dump_cache() const
{
    for (typename Cache<Vertex, GFRetVal>::const_iterator it = cache.begin(); it != cache.end(); ++it)
    {
        const Vertex& v1(it->first.v1);
        const Vertex& v2(it->first.v2);
        std::cout<<v1<<std::endl;
        std::cout<<v2<<std::endl;
    }
    return;
}

template <class Configuration, class GreensFunction, typename FPType>
void ObservableContainer<Configuration, GreensFunction, FPType>::measure(const Configuration& conf)
{
    Wick<Configuration, GreensFunction, FPType> wick;
    cache.template measure<Wick<Configuration, GreensFunction, FPType> > (wick, conf);
    DoWick<Vertex, GFRetVal> wickobj(cache);
    //now measure the actual observables
//the next pragma specifies that we enter a parallel region with the potential of #cores threads
#pragma omp parallel
//now we restrict that only one thread generates all threads
#pragma omp single
    //create the tasks for measuring
    for (int k = 0; k < static_cast<int>(container.size()); ++k)
    {
#ifdef _OPENMP
    std::cout<<"Generating task for measurement"<<std::endl;
//here is the actual work item of a thread
    #pragma omp task shared(conf, wickobj) untied
    {
    std::cout<<"measuring"<<std::endl;
#endif
        container[k]->evaluate(conf, wickobj);
#ifdef _OPENMP     
        std::cout<<"finished measurement!"<<std::endl;
    }
#endif
    }
    #pragma omp taskwait
    std::cout<<"exiting measure"<<std::endl;
}

template <class Configuration, class GreensFunction, typename FPType>
void ObservableContainer<Configuration, GreensFunction, FPType>::dryrun()
{
    //Note: via dryobj is the access to the cache
    DryRun<Vertex, GFRetVal> dryobj(cache);
    for (unsigned int k = 0; k < container.size(); ++k)
        container[k]->dryrun(dryobj);
    cache.init();
    return;
}

template <class Configuration, class GreensFunction, typename FPType>
ObservableBase<Configuration, typename GreensFunction::FreeGreensFunctionReturnValueType>*&
ObservableContainer<Configuration, GreensFunction, FPType>::operator[](unsigned int k) const
{
    return container[k];
}

template <class Configuration, class GreensFunction, typename FPType>
void ObservableContainer<Configuration, GreensFunction, FPType>::add_Cache_Value(Vertex& v1, Vertex& v2)
{
    cache.push_back(VertexPair<Vertex>(v1, v2));
}

template <class Configuration, class GreensFunction, typename FPType>
void ObservableContainer<Configuration, GreensFunction, FPType>::add_Observable(ObservableBase<Configuration, GFRetVal>* arg)
{
    container.push_back(arg);
}
#endif
