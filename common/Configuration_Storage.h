#ifndef CONFIGURATION_STORAGE_H
#define CONFIGURATION_STORAGE_H
#include "spin_symmetry_trait.h"

/**
This class is for determining the phase of a weight as from the prefactors of the partition function.
We declare it here, but the actual implementation is in ddqmc.h
0 is the phase of an imaginary time model.
1 is the phase of a real time model
*/
template<typename GreensFunction, class Configuration, unsigned int>
struct PhaseEvaluator
{
    static inline void determinePhase(Configuration& );
};

/**
This is a wrapper around the basic container for storing the Configuration. It's mostly there for storing
the value of the Determinant. The type of the Determinant(real/complex) is accessible via DeterminantType.
It mimics the behaviour of the standard STL Containers.
*/
template <class C, typename DeterminantType, typename FPTypeT, class MatrixType, bool Spin_Symmetric_model = false>
class Configuration_Storage
{
public:
    enum
    {
        Has_Spin_Symmetry = Spin_Symmetric_model
    };
    typedef C Container;///< the type of the used Container
    typedef DeterminantType DetType;///< the type of the used determinant
    typedef Matrix_Container<Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model> > MatConf;///< A typedef for the type of the StorageContainer for the Matrices
    typedef MatrixType MatType;///< the Type of the used matrices
    typedef typename Container::value_type value_type;///< the typedef for the values of the container
    typedef typename Container::const_iterator const_iterator;///< A const iterator
    typedef typename Container::const_reference const_reference;///< the type of the used reference
    typedef typename Container::iterator iterator;
    typedef typename Container::size_type size_type;
    typedef FPTypeT FPType;
    MatConf matcont;
    DeterminantType incphase;///< a variable for storing the incrementally changed phase
    DeterminantType phase;///< this is the phase of the Configuration. if the phase is real, this boils down to the sign of the Configuration
    inline Configuration_Storage(DetType init) : matcont(0), incphase(init), phase(init), initializer(init) {}
    inline Configuration_Storage() : incphase(1.0), phase(1.0), initializer(1.0) {}
    /**
    Copy constructor. YOU have to make sure that the Configuration_storage is in a valid state for copying!!!
    @param rhs the object to be copied
    */
    inline Configuration_Storage(const Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>& rhs) : matcont(rhs.matcont), incphase(rhs.incphase), phase(rhs.phase), container(rhs.container), initializer(rhs.initializer) {}
    /**
    A function to reset the scalar values to their initial values. Note that the matrices are not reset
    */
    inline void lazyreset() throw();
    /**
    A function to reset the matrices to some initial Matrices
    */
    inline void resetmatrices() throw();
    /**
    A function to copy some values of another container
    @param rhs the other Configuration_storage whose scalar values we copy. Note that the matrices are not copied over!
    */
    inline void lazyclone(const Configuration_Storage<C, DeterminantType, FPType, MatrixType, Spin_Symmetric_model>& rhs);
    /**
    return the number of Vertices stored in the container
    @return the size of the container
    */
    inline typename Container::size_type size() const
    {
        return container.size();
    }
    /**
    an accessor to the elements of the container
    @param k the Element that you want to access
    */
    inline value_type& operator[](unsigned int k)
    {
        return container[k];
    }
    /**
    an accessor to the elements of the container. the const version
    @param k the Element that you want to access
    */
    inline const value_type& operator[](unsigned int k) const
    {
        return container[k];
    }
    /**
    A function for obtaining a const_iterator to the beginning of the container
    @return a const_iterator to the beginning of the container
    */
    inline const_iterator begin() const
    {
        return container.begin();
    }
    inline iterator begin()
    {
        return container.begin();
    }
    inline iterator end()
    {
        return container.end();
    }
    inline const_iterator end() const
    {
        return container.end();
    }
    /**
    Store the argument in the container
    @param arg the element you want to append to the container
    */
    inline void push_back(const value_type& arg)
    {
        container.push_back(arg);
    }
    inline void pop_back()
    {
        container.pop_back();
    }
    inline iterator erase(iterator pos)
    {
        return container.erase(pos);
    }
    inline const_reference back() const
    {
        return container.back();
    }
    template<class Greensfunction>
    bool needsreinit(const Alpha<FPType>& alpha);
    /**
    A function to determine the product of the determinants. Useful for debugging
    @param alpha the alpha Object to use
    */
    template<class Greensfunction>
    DeterminantType dets(const Alpha<FPType>& alpha) const;
    template <class GreensFunction>
    void reevaluate(const Alpha<FPType>& alpha);
    template <class GreensFunction>
    inline typename GreensFunction::FreeGreensFunctionReturnValueType measure(const value_type& creator, const value_type& destructor) const;
//private:
    Container container;///< the Container in which the Vertices are stored
    DeterminantType initializer;///< the value we reset everything to
};

template <int Realtime>
class Keldysh_Brancher;

template <>
struct Keldysh_Brancher<0>
{
    template <class GreensFunction, class Mat, class Cont>
    static inline typename GreensFunction::FreeGreensFunctionReturnValueType measure(const Mat& mat, const Cont& cont, const typename GreensFunction::Vertex& creator, const typename GreensFunction::Vertex& destructor)
    {
      return mat.template measure<GreensFunction>(cont, creator, destructor);
    }
};

template <>
struct Keldysh_Brancher<1>
{
    template <class GreensFunction, class Mat, class Cont>
    static inline typename GreensFunction::FreeGreensFunctionReturnValueType measure(const Mat& mat, const Cont& cont, const typename GreensFunction::Vertex& creator, const typename GreensFunction::Vertex& destructor)
    {
        const typename GreensFunction::Vertex creatorc = GreensFunction::switchBranch(creator);
        const typename GreensFunction::Vertex destructorc = GreensFunction::switchBranch(destructor);
        return static_cast<typename GreensFunction::FreeGreensFunctionReturnValueType>(0.5)*(
                   mat.template measure<GreensFunction>(cont, creator, destructor) + mat.template measure<GreensFunction>(cont, creatorc, destructorc));
    }
};

template<class C, typename DeterminantType, typename FPTypeT, class MatrixType, bool Spin_Symmetric_model>
template<class GreensFunction>
typename GreensFunction::FreeGreensFunctionReturnValueType Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>::measure(const value_type& creator, const value_type& destructor) const
{
    return Keldysh_Brancher<GreensFunction::timeevolution>::template measure<GreensFunction, MatConf, Container>(matcont, container, creator, destructor);
}

template<class C, typename DeterminantType, typename FPTypeT, class MatrixType, bool Spin_Symmetric_model>
template<class GreensFunction>
void Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>::reevaluate(const Alpha<FPType>& alpha)
{
    matcont.template reevaluate<GreensFunction, Alpha<FPType> >(container, alpha, incphase);
    PhaseEvaluator<GreensFunction, Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>, GreensFunction::timeevolution>::determinePhase(*this);
    return;
}

template<class C, typename DeterminantType, typename FPTypeT, class MatrixType, bool Spin_Symmetric_model>
template<class GreensFunction>
bool Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>::needsreinit(const Alpha<FPType>& alpha)
{
    PhaseEvaluator<GreensFunction, Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>, GreensFunction::timeevolution>::determinePhase(*this);
    if (this->size() == 0) return false;//empty configuration, nothing to do
    return matcont.template needsreinit<GreensFunction, Alpha<FPType> >(container, alpha);
}

template<class C, typename DeterminantType, typename FPTypeT, class MatrixType, bool Spin_Symmetric_model>
void Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>::lazyclone(const Configuration_Storage<C, DeterminantType, FPType, MatrixType, Spin_Symmetric_model>& rhs)
{
    initializer = rhs.initializer;
    incphase = rhs.incphase;
    phase = rhs.phase;
    container = rhs.container;
    return;
}

template<class C, typename DeterminantType, typename FPTypeT, class MatrixType, bool Spin_Symmetric_model>
void Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>::resetmatrices() throw()
{
    this->matcont.reset();
}

template<class C, typename DeterminantType, typename FPTypeT, class MatrixType, bool Spin_Symmetric_model>
void Configuration_Storage<C, DeterminantType, FPTypeT, MatrixType, Spin_Symmetric_model>::lazyreset() throw()
{
    this->phase = this->initializer;
    return;
}

/**
A stream Operator for outputting the Configuration_Storage container
*/
template < template <class, class, class, class, bool> class Conf, template <class tp, class alloc = std::allocator<tp> > class Cont, class Out_Vertex, class det, class FPType, class MatrixType, bool Spin_Symmetry>
std::ostream& operator<<(std::ostream& out, const Conf<Cont<Out_Vertex>, det, FPType, MatrixType, Spin_Symmetry>& rhs )
{
    out<<"C_"<<rhs.size()<<" = [ "<<std::endl;
    typename Cont<Out_Vertex>::const_iterator it;
    for ( it = rhs.begin(); it != rhs.end(); it++)
        out<<(*it)<<", "<<std::endl;
    out<<"]";
    return out;
}
#endif
