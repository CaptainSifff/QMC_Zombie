/***************************************************************************
 *   *
 *  fgoth@physik.uni-wuerzburg.de *
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
#ifndef MEMARRAY_1D
#define MEMARRAY_1D
#include "meta.h"
#include "TypePromoter.h"
#include "MTL_Common.h"

namespace MTLICCL
{

template <
typename T_Elementtype>
struct Config1D
{
    typedef T_Elementtype ElementType;
    enum {
        Is_Static = false,
        Is_Dynamic = true,
        Is_OneD = true,
        rows = 0,
        columns = 0,
    };
};

template <class, class, class>//some forward declaration
class Multiplication_specialization;

template < class Generator >
class MemArray1D
{
public:
    typedef Generator Config;
    typedef typename Config::ElementType Elementtype;

private:
    unsigned int _rows;
    unsigned int _columns;
    Elementtype* MTL_RESTRICT arr;
friend class Multiplication_specialization<MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > >, MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > >, MTLICCL::Matrix<MTLICCL::MemArray1D<MTLICCL::Config1D<double> > > >;
public:
    inline unsigned int  Rows() const throw()
    {
        return _rows;
    }
    inline unsigned int Columns() const throw()
    {
        return _columns;
    }
    inline MemArray1D(unsigned int r , unsigned int c) : _rows(r), _columns(c)
    {
#ifndef __USE_XOPEN2K
        arr = new Elementtype[r*c];
#else
        int ret = posix_memalign((void**)(&arr), 64, sizeof(Elementtype)*r*c);
        if(unlikely(ret != 0))
        {
          throw std::bad_alloc();
        }
#endif
    }
    inline MemArray1D() : _rows(0), _columns(0), arr(NULL) {}
    inline MemArray1D( const MemArray1D&MTL_RESTRICT arg ) : _rows(arg.Rows()), _columns(arg.Columns())
    {
      const uint r = arg.Rows();
      const uint c = arg.Columns();
#ifndef __USE_XOPEN2K
        arr = new Elementtype[r*c];
#else
        int ret = posix_memalign((void**)(&arr), 64, sizeof(Elementtype)*r*c);
        if(unlikely(ret != 0))
        {
          throw std::bad_alloc();
        }
#endif
        for ( unsigned int k = 0 ; k < r ; ++k )
            for ( unsigned int i = 0 ; i < c ; ++i )
                SetElement( k, i, arg(k,i));
        return;
    }
    inline MemArray1D& operator=(const MemArray1D&MTL_RESTRICT rhs)
    {//check for self-assignment removed, we assume the matrix class does that
	  const uint r = rhs.Rows();
	  const uint c = rhs.Columns();
            if ((this->Rows() != r) || (this->Columns() != c))
            {//If we copy around matrices of the same size we don't need to reallocate
                _rows = r;
                _columns = c;
#ifndef __USE_XOPEN2K
		 delete [] arr;
                 arr = new Elementtype[r*c];
#else
		 free(arr);
                 int ret = posix_memalign((void**)(&arr), 64, sizeof(Elementtype)*r*c);
                 if(unlikely(ret != 0))
                 {
                    throw std::bad_alloc();
                 }
#endif
            }
             for ( unsigned int k = 0 ; k < r; ++k )//likely optimizable by memcpy
                 for ( unsigned int i = 0 ; i < c; ++i )
                     SetElement( k, i, rhs(k,i));
        return *this;
    }
    template <typename T>
    inline MemArray1D& operator=(const Matrix<T>&MTL_RESTRICT rhs)
    {
    //No check for self-assignment since we are distinct types
            const unsigned int r = rhs.Rows();
            const unsigned int c = rhs.Columns();
            if ((r != this->Rows()) || (c != this->Columns()))
            {//If we copy around matrices of the same size we don't need to reallocate
                _rows = r;
                _columns = c;
#ifndef __USE_XOPEN2K
		 delete [] arr;
                 arr = new Elementtype[r*c];
#else
		 if(arr != NULL)
		   free(arr);
                 int ret = posix_memalign((void**)(&arr), 64, sizeof(Elementtype)*r*c);
                 if(unlikely(ret != 0))
                 {
                    throw std::bad_alloc();
                 }
#endif
            }
            for ( unsigned int k = 0 ; k < r; ++k )
                for ( unsigned int i = 0 ; i < c; ++i )
                    SetElement( k, i, rhs(k,i));
        return *this;
    }
#ifdef HAS_RVALUE_REFERENCES
    inline MemArray1D(MemArray1D&& arg ) : _rows(arg.Rows()), _columns(arg.Columns())
    {
        arr = arg.arr;
        arg.arr = NULL;
        return;
    }
    inline MemArray1D& operator=(MemArray1D&& rhs)
    {
        if (this != &rhs)
        {
            _rows = rhs.Rows();
            _columns = rhs.Columns();
#ifndef __USE_XOPEN2K
            delete [] arr;
#else
	    if(arr != NULL)
	      free(arr);
#endif
            arr = rhs.arr;
            rhs.arr = NULL;
        }
        return *this;
    }
#endif
    inline ~MemArray1D()
    {
#ifndef __USE_XOPEN2K
            delete [] arr;
#else
	    if(arr != NULL)
	      free(arr);
#endif
    }
    inline const Elementtype& operator() (unsigned int i , unsigned int k) const MTL_RESTRICT throw()
    {//const Version
        return arr[i * Columns() + k];
    }
    inline Elementtype& operator() (unsigned int i , unsigned int k) throw()
    {//Non-const Version
        return arr[i * Columns() + k];
    }
    inline const Elementtype& GetElement (unsigned int i , unsigned int k) const MTL_RESTRICT throw()
    {
        return arr[i * Columns() + k];
    }
    inline void SetElement (unsigned int i , unsigned int k, Elementtype val) throw()
    {
        arr[i * Columns() + k] = val;
        return;
    }
    /**
    A function for resizing the Matrix while keeping the old entries
    */
    void resize(unsigned int newsize);
};

template <class Config>
void MemArray1D<Config>::resize(unsigned int newsize)
{
    //fill up...
    return;
}
}
#endif
