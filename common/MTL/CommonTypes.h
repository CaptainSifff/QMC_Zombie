/***************************************************************************
 *   Copyright (C) 2007-2013 by Florian Goth   *
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
//Header to define shortcuts for various common Matrix Types
#include <complex>

typedef MTLICCL::Matrix< MTLICCL::Static<MTLICCL::Config<float,4,4> > > Matrix4x4;
typedef MTLICCL::Matrix< MTLICCL::Static<MTLICCL::Config<float,3,3> > > Matrix3x3;
typedef MTLICCL::Matrix< MTLICCL::Static<MTLICCL::Config<float,4,1> > > Vec4;
//typedef MTLICCL::Matrix< MTLICCL::Static<MTLICCL::Config<float,3,1> > > Vec3;
typedef MTLICCL::Matrix< MTLICCL::Dynamic< MTLICCL::Config_Dynamic<double> > > Mat_Dynamic;

typedef MTLICCL::Matrix< MTLICCL::Dynamic< MTLICCL::Config_Dynamic<std::complex<double> > > > CmplxMat;
//typedef MTLICCL::Matrix< MTLICCL::MemArray1D<MTLICCL::Config1D<std::complex<double> > > > CmplxMat;