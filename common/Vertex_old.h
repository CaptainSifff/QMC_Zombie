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
#ifndef VERTEX_H
#define VERTEX_H
#include <iostream>
#include "MTL/common.h"

class Basic_Vertex
{
public:
    inline Basic_Vertex(float t1, int s) throw();
    float tau;
    int spin;
    inline bool operator<(const Basic_Vertex& rhs) const throw();
    inline bool operator==(const Basic_Vertex& rhs) const throw();
private:
};

Basic_Vertex::Basic_Vertex(float t1, int s) throw() : tau(t1), spin(s) {}

bool Basic_Vertex::operator==(const Basic_Vertex& rhs) const throw()
{
    return (std::fabs(tau - rhs.tau) < tiny<float>() ) && (spin == rhs.spin);
}

bool Basic_Vertex::operator<(const Basic_Vertex& rhs) const throw()
{
    if (tau >= rhs.tau)//if tau is larger... we return true, else we enter the if-body
    {
        if (std::fabs(tau - rhs.tau) < tiny<float>())//if the taus are equal we look at the spin
        {
            return spin < rhs.spin;
        }
        return false;
    }
    return true;
}

/**
The Basic Vertex of the Hubbard-Model
*/
class Hubbard_Vertex
{
public:
    /**
    Constructor of the Vertex.
    @param t1 the first imaginary time.
    @param s1 the first lattice site
    @param s the spin
    */
    inline Hubbard_Vertex(float t1, int s1, int s) throw();
    float tau;///< the imaginary time of the Vertex
    int site;///< the site of the Vertex
    int spin;///< the spin of the Vertex
    /**
    An operator that compares two vertices lexicographically.
    @param rhs the other vertex to compare against.
    @return true if this vertex is smaller than the other else false
    */
    inline bool operator<(const Hubbard_Vertex& rhs) const throw();
    /**
    a operator that probes the identity of two Vertices
    @param rhs the Vertex to compare against
    @return true if the Vertices are equal in all of its members, else false.
    */
    inline bool operator==(const Hubbard_Vertex& rhs) const throw();
private:
};

bool Hubbard_Vertex::operator==(const Hubbard_Vertex& rhs) const throw()
{
    return (std::fabs(tau - rhs.tau) < tiny<float>() ) && (site == rhs.site) && (spin == rhs.spin);
}

bool Hubbard_Vertex::operator<(const Hubbard_Vertex& rhs) const throw()
{
    if (tau >= rhs.tau)
    {
        if (std::fabs(tau - rhs.tau) < tiny<float>())
        {
            if (site >= rhs.site)
            {
                if (site == rhs.site)
                {
                    return spin < rhs.spin;
                }
                return false;
            }
            return true;
        }
        return false;
    }
    return true;
}

/**
An overloaded stream operator for printing the basic vertex to stdout
*/
inline std::ostream& operator<<(std::ostream& out, const Hubbard_Vertex& rhs )
{
    out<<"("<<rhs.site<<", "<<rhs.tau<<", "<<rhs.spin<<")";
    return out;
}

Hubbard_Vertex::Hubbard_Vertex(float t, int si, int s) throw(): tau(t), site(si), spin(s)
{
}
#endif
