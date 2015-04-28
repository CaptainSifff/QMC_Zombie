/***************************************************************************
 *   Copyright (C) 2006 by Florian Goth   *
 *   CaptainSifff@gmx.de   *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#ifndef BOUNDSCHECKER_H
#define BOUNDSCHECKER_H

#include "Memory.h"

namespace MTLICCL
{
template < class Memory >
class BoundsChecker : public Memory
{
private:
protected:
    void checkbounds(unsigned int i, unsigned int k) const
    {
        if ( ( i >= Memory::Rows() ) || (k >= Memory::Columns() ) )//>= because of C-Style Element-addressing
            throw("Index Out of Bounds!!!" );
    }
public:
    inline BoundsChecker(unsigned int i , unsigned int k) : Memory(i,k)
{}
    inline BoundsChecker()
    {}
    typedef typename Memory::Config Config;
    typedef typename Config::ElementType Elementtype;

    inline const Elementtype& GetElement (unsigned int i , unsigned int k) const
    {
        checkbounds(i,k);
        return Memory::GetElement ( i , k );
    }
    inline void SetElement (unsigned int i , unsigned int k, Elementtype val)
    {
        checkbounds(i,k);
        return Memory::SetElement ( i , k, val );
    }
};
};
#endif
