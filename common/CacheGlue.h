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
#ifndef CACHE_GLUE_H
#define CACHE_GLUE_H
#include <map>
#include <utility>
#include <stdexcept>
#include <valarray>
#include "generalPhysics.h"
#include "toFPType.h"
#include "MTL/MTL_Macros.h"

/**
This structure encapsulates two Vertices and provides comparison operators
*/
template <class Vertex>
struct VertexPair
{
    Vertex v1;///< the first vertex
    Vertex v2;///< the second vertex
    /**
    Constructor for the vertex-pair
    @param a the first vertex
    @param b the second vertex
    */
    VertexPair(Vertex a, Vertex b) : v1(a), v2(b) {}
    VertexPair(const VertexPair& rhs) : v1(rhs.v1), v2(rhs.v2) {}
    /**
    A comparison operator that probes for identity
    @param arg the other vertexpair to compare against
    @return true if the vertices in each VertexPair are the identical, else false
    */
    inline bool operator==(const VertexPair& arg) const;
    /**
    A smaller-than operator to compare two vertices lexicographically
    @param arg the other Vertex-Pair to compare against
    @return true if the 'string' (v1,v2) is lexicographcally smaller than the same string arg. Else return false
    */
    inline bool operator<(const VertexPair& arg) const;
};

template <class T>
inline VertexPair<T> make_Vertex_pair(T a, T b)
{
    return VertexPair<T>(a,b);
}

template <class Configuration, class GreensFunction, typename FPType>
class Wick
{
public:
    typedef typename GreensFunction::FreeGreensFunctionReturnValueType RetType;
    typedef typename Configuration::value_type Vertex;
    typedef Configuration Config;
    typedef GreensFunction GF;
/**
This function forwards the measurement to the Matrix-Container. The Matrix-Container can better decide which optimizations to make.
*/
    inline RetType measure(const Configuration& configuration, const typename GreensFunction::Vertex& creator, const typename GreensFunction::Vertex& destructor) const
    {
      RetType retval = configuration.template measure<GreensFunction>(creator, destructor);
      return retval;
    }
};

template <class Vertex, class GFRetVal>
class DryRun;//forward declare DryRun

template <class Vertex, class GFRetVal>
class Cache
{
public:
    typedef typename std::map<VertexPair<Vertex>, GFRetVal >::iterator iterator;
    typedef typename std::map<VertexPair<Vertex>, GFRetVal >::const_iterator const_iterator;
    typedef std::valarray<std::pair<const VertexPair<Vertex>, GFRetVal>* > LinearType;
    inline std::size_t size() const throw()
    {
        return cache.size();
    }
    inline GFRetVal get_key(const VertexPair<Vertex> vp) const MTL_PURE_FUNCTION
    {
      const_iterator it = cache.find(vp);
      if(unlikely(it == cache.end()))
      {
	std::cout<<"Vertices not found: "<<std::endl;
	std::cout<<vp.v1<<std::endl;
	std::cout<<vp.v2<<std::endl;
	std::cout<<"Aborting... Do something!!!!!"<<std::endl;
	exit(-1);
      }
        return it->second;
    }
    /**
    A function for obtaining a const_iterator to the beginning of the map
    @return a const_iterator to the beginning of the map
    */
    inline const_iterator begin() const
    {
        return cache.begin();
    }
    inline iterator begin()
    {
        return cache.begin();
    }
    inline iterator end()
    {
        return cache.end();
    }
    inline const_iterator end() const
    {
        return cache.end();
    }
    inline LinearType linear_array()
    {
        LinearType ret(this->size());
        int u = 0;
        for (typename std::map<VertexPair<Vertex>, GFRetVal>::iterator it = cache.begin(); it != cache.end(); ++it, ++u)
            ret[u] = &(*it);
        return ret;
    }
    inline void init()
    {
#ifdef _OPENMP
        items.resize(this->size());
        int u = 0;
        for (typename std::map<VertexPair<Vertex>, GFRetVal>::iterator it = cache.begin(); it != cache.end(); ++it, ++u)
            items[u] = &(*it);
#endif
    }
    template <class T>
    inline void measure(const T& wick, const typename T::Config& conf)
    {
#ifdef _OPENMP
#pragma omp parallel for
        for (int u = 0; u < items.size(); ++u)
        {
            const Vertex& v1(items[u]->first.v1);
            const Vertex& v2(items[u]->first.v2);
            items[u]->second = wick.measure(conf, v1, v2);
        }
#else
        for (typename std::map<VertexPair<Vertex>, GFRetVal >::iterator it = cache.begin(); it != cache.end(); ++it)
        {
            const Vertex& v1(it->first.v1);
            const Vertex& v2(it->first.v2);
            it->second = wick.measure(conf, v1, v2);
        }
#endif
    }
    inline Cache(unsigned int) {}
private:
    friend class DryRun<Vertex, GFRetVal>;
    inline void insert(VertexPair<Vertex>& arg)
    {
        cache.insert(std::make_pair(arg, GFRetVal(0.0)));
    }
    std::map<VertexPair<Vertex>, GFRetVal> cache;///< the cache with all the requested greensfunction values
#ifdef _OPENMP
    LinearType items;///< some temporary space where we store the pointers into the map in a linear fashion
#endif
};

template <class GFRetVal, typename FPType>
class Cache<Hubbard_Vertex<FPType>, GFRetVal>;

template <class GFRetVal, typename FPType>
class Cache_Iterator : std::iterator<std::input_iterator_tag, std::pair<VertexPair<Basic_Vertex<FPType> >, GFRetVal> >
{
public:
    inline Cache_Iterator& operator++()// equals ++it
    {
        ++it;
        if (it == cache[idx].end())
        {
            ++idx;
            it = cache[idx].begin();
        }
        return *this;
    }
    inline bool operator==(const Cache_Iterator& rhs) const
    {
        return (cache == rhs.cache) && (idx == rhs.idx) && (it == rhs.it);
    }
    inline bool operator!=(const Cache_Iterator& rhs) const
    {
        return !this->operator==(rhs);
    }
    inline std::pair<VertexPair<Basic_Vertex<FPType> >, GFRetVal>& operator*()
    {
        return *it;
    }
private:
    unsigned int idx;
    typename std::map<VertexPair<Basic_Vertex<FPType> >, GFRetVal >::iterator it;
    Cache<Hubbard_Vertex<FPType>, GFRetVal>& cache;
    Cache_Iterator(Cache<Hubbard_Vertex<FPType>, GFRetVal>& arg) : idx(0), cache(arg), it( cache[idx].begin() ) {}
};

template <class GFRetVal, typename FPType>
class Cache<Hubbard_Vertex<FPType>, GFRetVal>
{
public:
    typedef Cache_Iterator<GFRetVal, FPType> iterator;
    typedef const Cache_Iterator<GFRetVal, FPType> const_iterator;

    inline std::size_t size() const throw() MTL_PURE_FUNCTION
    {
        std::size_t ret = 0;
        for (unsigned int k = 0; k < chain_len * chain_len; ++k)
            ret+= cache[k].size();
        return ret;
    }
    inline GFRetVal get_key(const VertexPair<Hubbard_Vertex<FPType> > vp) const
    {
        return access(vp.v1.site, vp.v2.site).find(make_Vertex_pair(Basic_Vertex<FPType>(vp.v1.tau, vp.v1.spin), Basic_Vertex<FPType>(vp.v2.tau, vp.v2.spin)))->second;
    }
    inline const_iterator begin() const
    {
        return cache.begin();
    }
    inline iterator begin()
    {
        return cache.begin();
    }
    inline iterator end()
    {
        return cache.end();
    }
    inline const_iterator end() const
    {
        return cache.end();
    }
    inline Cache(const unsigned int len) : cache(new std::map<VertexPair<Basic_Vertex<FPType> >, GFRetVal>[len * len]),
            chain_len(len)
    {}
    inline ~Cache()
    {
        delete [] cache;
    }
    /**
    A function to measure the greensfunction values that are stored in the cache
    @param wick The object that performs the measurement on a Configuration
    @param conf The configuration
    */
    template <class W>
    inline void measure(const W& wick, const typename W::Config& conf)
    {
            std::cout<<"Beginning measurement of Greensfunctions"<<std::endl;
//somehow playing with the chunk size could probably help with the scheduling...
#pragma omp parallel for schedule(dynamic)
        for (unsigned int k = 0; k < chain_len*chain_len; ++k)
        {
//            for (unsigned int j = 0; j < chain_len; ++j)
//            {
	  const uint idx_v1 = k/chain_len;
	  const uint idx_v2 = k%chain_len;
          std::map<VertexPair<Basic_Vertex<FPType> >, GFRetVal>& lm(access(idx_v1, idx_v2));
                for (typename std::map<VertexPair<Basic_Vertex<FPType> >, GFRetVal>::iterator it = lm.begin(); it != lm.end(); ++it)
                {
                    const VertexPair<Basic_Vertex<FPType> >& vp(it->first);
                    const Hubbard_Vertex<FPType> v1(idx_v1, vp.v1.tau, vp.v1.spin);
                    const Hubbard_Vertex<FPType> v2(idx_v2, vp.v2.tau, vp.v2.spin);
                    it->second = wick.measure(conf, v1, v2);
                }
//            }
        }
        std::cout<<"measurement of Greensfunctions done!"<<std::endl;
    }
    inline void init()
    {
    }
private:
    inline void insert(VertexPair<Hubbard_Vertex<FPType> >& arg)
    {
        GFRetVal zero = 0;
        access(arg.v1.site, arg.v2.site).insert(std::make_pair(make_Vertex_pair(Basic_Vertex<FPType>(arg.v1.tau, arg.v1.spin), Basic_Vertex<FPType>(arg.v2.tau, arg.v2.spin)), zero));
    }
    friend class DryRun<Hubbard_Vertex<FPType>, GFRetVal>;
    std::map<VertexPair<Basic_Vertex<FPType> >, GFRetVal>* cache;///< the cache with all the requested greensfunction values
    unsigned int chain_len;
    inline const std::map<VertexPair<Basic_Vertex<FPType> >, GFRetVal>& access(unsigned int i, unsigned int k) const
    {
        return cache[i*chain_len + k];
    }
    inline std::map<VertexPair<Basic_Vertex<FPType> >, GFRetVal>& access(unsigned int i, unsigned int k)
    {
        return cache[i*chain_len + k];
    }
};

/**
We need to provide the Greensfunctionvalue-cache with the necessary data which Greensfunction values we want to evaluate at every Configuration. For this, each observable provides a dryrun() method that takes a DryRun object which takes care of storing the Objects into the cache.
*/
template <class Vertex, class GFRetVal>
class DryRun
{
public:
    typedef typename SignToFPType<GFRetVal>::Ret FPType;
    /**
    Constructor of the DryRun Object
    @param arg a reference to the storage of the cache
    */
    inline DryRun(Cache<Vertex, GFRetVal>& arg) : cache(arg) {}
    /**
    This ()-operator takes the two Vertices at which to evaluate the Greensfunction and stores them in the cache
    @param v1 the first Vertex (the adjungated fermi-operator)
    @param v2 the second Vertex
    */
    inline void operator()(Vertex v1, Vertex v2);
    /**
    A function that maps the previously used doWickonSector to the Cache architecture
    the template parameter selects the Spin sector
    @param site_i the site of the first Vertex
    @param tau_i the time of the first vertex
    @param site_j the site of the second vertex
    @param tau_j the time of the second vertex
    */
    template<SPINS Spin>
    inline void onSector(const int site_i, const FPType tau_i, const int site_j, const FPType tau_j);
private:
    //declare certain things private to not allow copying this object
    DryRun();
    DryRun(const DryRun&);
    DryRun& operator=(const DryRun&);
    Cache<Vertex, GFRetVal>& cache;///< A reference to the structure for the cache
};

/**
This class encapsulates the access from an observable to the cache. The Observables provide evaluate() methods that take this object, which will be used for evaluating the Greensfunctionvalues.
*/
template <class Vertex, class GFRetVal>
class DoWick
{
public:
    typedef typename SignToFPType<GFRetVal>::Ret FPType;
    /**
    The constructor of the DoWick object
    @param arg a reference to the object that has the cached values
    */
    inline DoWick(const Cache<Vertex, GFRetVal>& arg) : cache(arg) {}
    /**
    This ()-operator takes the two Vertices at which to evaluate the Greensfunction and stores them in the cache
    @param v1 the first Vertex (the adjungated fermi-operator)
    @param v2 the second Vertex
    */
    inline GFRetVal operator()(Vertex v1, Vertex v2) const;
    /**
    This function mimicks the behaviour of the previously used doWickonSector function. It returns the value of the Greensfunction for the given values
    @param site_i the site of the first Vertex
    @param tau_i the time of the first vertex
    @param site_j the site of the second vertex
    @param tau_j the time of the second vertex
    @return the value of the Greensfunction for the given values
    */
    template<SPINS Spin>
    inline GFRetVal onSector(const int site_i, const FPType tau_i, const int site_j, const FPType tau_j) const;
private:
    //declare certain things private to not allow copying this object
    DoWick();
    DoWick(const DoWick&);
    DoWick& operator=(const DoWick&);
    const Cache<Vertex, GFRetVal>& cache;
};

template <class Vertex, typename FPType>
struct DryRunHelper
{
    static inline VertexPair<Vertex> createPair(const int site_i, const FPType tau_i, const int site_j, const FPType tau_j, SPINS spin);
};

template <typename FPType>
struct DryRunHelper<Basic_Vertex<FPType>, FPType>
{
    static inline VertexPair<Basic_Vertex<FPType> > createPair(const int site_i, const FPType tau_i, const int site_j, const FPType tau_j, SPINS spin)
    {
        return VertexPair<Basic_Vertex<FPType> >(Basic_Vertex<FPType>(tau_i, spin), Basic_Vertex<FPType>(tau_j, spin));
    }
};

template <typename FPType>
struct DryRunHelper<Hubbard_Vertex<FPType>, FPType>
{
    static inline VertexPair<Hubbard_Vertex<FPType> > createPair(const int site_i, const FPType tau_i, const int site_j, const FPType tau_j, SPINS spin)
    {
        return VertexPair<Hubbard_Vertex<FPType> >(Hubbard_Vertex<FPType>(site_i, tau_i, spin), Hubbard_Vertex<FPType>(site_j, tau_j, spin));
    }
};

template< class Vertex, class GFRetVal>
template<SPINS Spin>
GFRetVal DoWick<Vertex, GFRetVal>::onSector(const int site_i, const FPType tau_i, const int site_j, const FPType tau_j) const
{
    return cache.get_key(DryRunHelper<Vertex, FPType>::createPair(site_i, tau_i, site_j, tau_j, Spin));
}

template< class Vertex, class GFRetVal>
template<SPINS Spin>
void DryRun<Vertex, GFRetVal>::onSector(const int site_i, const FPType tau_i, const int site_j, const FPType tau_j)
{
    VertexPair<Vertex> ce = DryRunHelper<Vertex, FPType>::createPair(site_i, tau_i, site_j, tau_j, Spin);
    this->cache.insert(ce);
}

template< class Vertex, class GFRetVal>
GFRetVal DoWick<Vertex, GFRetVal>::operator()(Vertex v1, Vertex v2) const
{
    VertexPair<Vertex> ce(v1, v2);
    return cache.get_key(ce);
}

template< class Vertex, class GFRetVal>
void DryRun<Vertex, GFRetVal>::operator()(Vertex v1, Vertex v2)
{
    VertexPair<Vertex> ce(v1, v2);
    this->cache.insert(ce);
}

template <class Vertex>
bool VertexPair<Vertex>::operator<(const VertexPair& arg) const
{
    if (!(v1 < arg.v1))
    {
        //greater-than-branch
        if (v1 == arg.v1)
        {
            return v2 < arg.v2;
        }
        return false;
    }
    return true;
}

template <class Vertex>
bool VertexPair<Vertex>::operator==(const VertexPair& arg) const
{
    return (arg.v1 == v1) && (arg.v2 == v2);
}
#endif
