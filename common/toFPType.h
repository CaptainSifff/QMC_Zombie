#ifndef TOFPTYPE_H
#define TOFPTYPE_H
#include <complex>
template <typename T> struct SignToFPType
{
    inline T zero() const throw() {
        return 0;
    }
    typedef T Ret;
};

template <typename S>
struct SignToFPType<std::complex<S> >
{
    inline std::complex<S> zero() const throw() {
        return std::complex<S>(0.0, 0.0);
    }
    typedef S Ret;
};
#endif
