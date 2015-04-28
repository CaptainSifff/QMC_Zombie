/***************************************************************************
 *   Copyright (C) 2013 by Florian Goth   *
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
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <vector>

#include "Matrix.h"
#include "MatrixGenerator.h"
#include "CommonTypes.h"
#include "MemArray1D.h"

using namespace std;

/*Small tool to pass a text file containing a matrix to our eigenvalue solver*/
int main(int argc, char *argv[])
{
  typedef double FPType;
  std::ifstream fi(argv[1]);
  std::vector<std::vector<FPType> > matinput;
  std::string line;
  while(getline(fi, line, '\n'))
  {
    matinput.push_back(std::vector<FPType>());
    istringstream iss(line);//iss now contains the line of strings
    double d;
    while (iss >> d) matinput.back().push_back(d);
  }
Mat_Dynamic mat(matinput.size(), matinput.size());
for(uint i = 0; i < matinput.size(); ++i)
  for(uint k = 0; k < matinput.size(); ++k)
    mat(i, k) = matinput[i][k];
MTLICCL::EigenSystemSolver<Mat_Dynamic> ess(mat);
ess.calculateEVs();
Mat_Dynamic z(ess.tridiagonalize());
std::vector<MTLICCL::SolutionPair<double> > esystem = ess.eigensystem();
for(uint k = 0; k < esystem.size(); ++k)
  esystem[k].print();
// Mat_Dynamic evscplx(dim, dim);
//  for(uint k = 0; k < dim; ++k)
//    for(uint j = 0; j < dim; ++j)
//    {
//      evscplx(j, k) = esystem[k].evector[j];
//    }
//    std::cout<<z*evscplx<<std::endl;

/*
CmplxMat testmat(dim, dim);
for(uint k = 0; k < dim; ++k)
{
  testmat(k,k) = 2.0;
  if(k+1<dim)
  {
  testmat(k,k+1) = 1.0;
  testmat(k+1,k) = 1.0;
  }
  if(k+2 < dim)
  {
  testmat(k,k+2) = std::complex<double>(0.0, 1.0)*0.5;
  testmat(k+2,k) = 0.5*std::complex<double>(0.0, -1.0);
  }
  if(k+3 < dim)
  {
  testmat(k,k+3) = std::complex<double>(1.0,+ 1.0);
  testmat(k+3,k) = std::complex<double>(1.0, -1.0);
  }
}
std::cout<<testmat<<std::endl;
MTLICCL::EigenSystemSolver<CmplxMat> ess(testmat);
MTLICCL::EigenSystemSolver<CmplxMat> ess2(testmat);
ess2.calculateEVs();
CmplxMat zcplx(ess2.tridiagonalize());
//CmplxMat zcplx(ess.tridiagonalize());
std::vector<MTLICCL::SolutionPair<double> > esystemcplx = ess2.eigensystem();
 CmplxMat evscplx(zcplx.Rows(), zcplx.Rows());
 for(uint k = 0; k < zcplx.Rows(); ++k)
   for(uint j = 0; j < zcplx.Rows(); ++j)
   {
     evscplx(j, k) = esystemcplx[k].evector[j];
   }
   std::cout<<zcplx*evscplx<<std::endl;
 CmplxMat zhc(~(zcplx*evscplx));
 for(uint k = 0; k < zhc.Rows(); ++k)
   for(uint i = 0; i < zhc.Columns(); ++i)
     zhc(k,i) = conj(zhc(k,i));
   CmplxMat id(zhc*(zcplx*evscplx));
 for(uint k = 0; k < id.Rows(); ++k)
 {
     double sum = 0.0;
   for(uint i = 0; i < id.Rows(); ++i){
     if(i!=k)
     sum += abs(id(i,k));
   }
   std::cout<<"Zeilensumme "<<k<<"    : "<<scientific<<sum<<" + "<<abs(id(k,k))<<std::endl;
 }
 std::cout.precision(12);
 ofstream file("evs1000.txt");
 file.precision(12);
 for(uint k = 0; k < esystemcplx.size(); ++k)
   file<<esystemcplx[k].evalue<<", ";*/
//   for(uint j = 0; j < esystemcplx.size(); ++j)
//     file<<k<<" "<<esystemcplx[j].evalue<<std::endl;
//   
//  CmplxMat evscplx(zcplx.Rows(), zcplx.Rows());
//  for(uint i = 0; i < zcplx.Rows(); ++i)
//    for(uint j = 0; j < zcplx.Rows(); ++j)
//    {
//      evscplx(j, i) = esystemcplx[i].evector[j];
//    }
// //   std::cout<<zcplx*evscplx<<std::endl;
//  CmplxMat zhc(~(zcplx*evscplx));
//  for(uint j = 0; j < zhc.Rows(); ++j)
//    for(uint i = 0; i < zhc.Columns(); ++i)
//      zhc(j,i) = conj(zhc(j,i));
//    CmplxMat id(zhc*(zcplx*evscplx));
//  std::cout<<id<<std::endl;
//  for(uint j = 0; j < id.Rows(); ++j)
//  {
//      double sum = 0.0;
//    for(uint i = 0; i < id.Rows(); ++i){
//      if(i!=j)
//      sum += abs(id(i,j));
//    }
//     std::cout.precision(12);
//    summen<<"Zeilensumme "<<j<<"    : "<<scientific<<sum<<" + "<<abs(id(k,k))<<std::endl;
//  }
    return EXIT_SUCCESS;
}
