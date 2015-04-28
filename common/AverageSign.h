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
#ifndef AVERAGESIGN_H
#define AVERAGESIGN_H
#include "Zombie.h"
#include "scalar_observable.h"
/**
A class for measuring the average Sign.
*/
template <class Configuration, typename SignType>
class AverageSign : public mc_analysis::Scalar_Observable<SignType>
{
public:
    /**
    The Constructor for the Average Sign
    */
    inline AverageSign()
    {}
    inline ~AverageSign(){}
    /**
    This determines the Average sign for the given configuration
    @param configuration the configuration
    */
    inline void evaluate(const Configuration&) throw();
    template <class Net>
    inline void sendData(Net&);
private:
};

template <class Configuration, typename SignType>
template <class Net>
void AverageSign<Configuration, SignType>::sendData(Net& z)
{
    this->rebin_data();
    z.template sendtoserver<SignType>(this->mean());
}

template <class Configuration, typename SignType>
void AverageSign<Configuration, SignType>::evaluate(const Configuration& configuration) throw()
{
    this->push_back(configuration.phase);
    return;
}
#endif
