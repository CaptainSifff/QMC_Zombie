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
#ifndef MPICLIENT_H
#define MPICLIENT_H
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#include <stdint.h>
#include <valarray>
#include <complex>
#include <stdexcept>
#include "MPITypeDeduction.h"


class MPIClient
{
public:
    /**
    Constructor call. it needs argc and argv from the command line, to pass it on to MPI
    @param argc argc as given from the main funtion
    @param argv argv as given from the main function
    */
    inline MPIClient(int argc, char *argv[]);
    template <typename T>
    inline void sendtoserver(const T&);
    template <typename T>
    inline T recv();
    /**
    Blocking receive. This function tries to receive data. Note that it blocks all further execution.
    @return the requested value
    */
    template <typename T>
    inline T recvblock();
    /**
    This function checks whether there is data sent to us.
    @return true if there is incoming data
    */
    inline bool isdataready();
    inline ~MPIClient();
private:
    MPI_Status status;
    bool status_is_uptodate;
    void status_valid();
};

template<typename T>
T MPIClient::recvblock()
{
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    status_is_uptodate = true;
    return this->recv<T>();//in MPI recv is a blocking operation
}

void MPIClient::status_valid()
{
    if (!status_is_uptodate) throw(std::logic_error("[MPIClient] Status object is not valid!!!"));
}

bool MPIClient::isdataready()
{
    //This function immediately returns
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    status_is_uptodate = true;
    return flag;
}

template <>
inline std::string MPIClient::recv<std::string>()
{
    status_valid();
    int length;
    MPI_Get_count(&status, MPI_CHAR, &length);
    char* buf = new char[length];
    MPI_Recv(buf, length, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    status_is_uptodate = false;
    return std::string(buf, length);
}

template <>
inline std::valarray<char> MPIClient::recv<std::valarray<char> >()
{
    status_valid();
    int length;
    MPI_Get_count(&status, MPI_CHAR, &length);
    std::valarray<char> retval;
    retval.resize(length);
    MPI_Recv(&(retval[0]), length, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    status_is_uptodate = false;
    return retval;
}

template <typename T>
inline T MPIClient::recv()
{
    status_valid();
    T retval;
    status_valid();
    MPI_Recv(&retval, 1, mpi::MPIType<T>::get(), status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    status_is_uptodate = false;
    return retval;
}

template <typename T>
inline void MPIClient::sendtoserver(const T& x)
{
    T t(x);
    MPI_Send(&t, 1, mpi::MPIType<T>::get(), 0, 1, MPI_COMM_WORLD);
}

template <>
inline void MPIClient::sendtoserver<std::valarray<double> >(const std::valarray<double>& x)
{
    std::cout<<"sending to server"<<std::endl;
    std::valarray<double> t(x);
    MPI_Send(&t[0], x.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
}

template <>
inline void MPIClient::sendtoserver<std::valarray<std::complex<double> > >(const std::valarray<std::complex<double> >& x)
{
//this implementation uses the OpenMPI MPI_DOUBLE_COMPLEX define
    std::cout<<"sending to server"<<std::endl;
    std::valarray<std::complex<double> > t(x);
    MPI_Send(&t[0], x.size(), MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
}

template <>
inline void MPIClient::sendtoserver<std::valarray<std::valarray<std::complex<double> > > >(const std::valarray<std::valarray<std::complex<double> > >& x)
{
//this implementation uses the OpenMPI MPI_DOUBLE_COMPLEX define
    uint32_t nrof_functions = x.size();
    MPI_Request request;
    if (MPI_Isend(&nrof_functions, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, &request) != 0)//send the number of the transmitted functions to the server
    {
        std::cout<<"An error happened while transmitting the number of functions!"<<std::endl;
        exit(-1);
    }
    std::cout<<"sending to server"<<std::endl;
    std::complex<double>* arr = new std::complex<double>[nrof_functions * x[0].size()];
    //copy into an internal buffer
    for (unsigned int k = 0; k < nrof_functions; ++k)
        for (std::size_t j = 0; j < x[0].size(); ++j)
            arr[k*x[0].size() + j] = x[k][j];
    //send the whole buffer at once
    MPI_Status st;
    MPI_Wait(&request, &st);
    if (MPI_Send(arr, nrof_functions * x[0].size(), MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD) != 0)
    {
        std::cout<<"An error happened while transmitting the actual data!"<<std::endl;
        exit(-1);
    }
    delete [] arr;
}

template <>
inline void MPIClient::sendtoserver<std::valarray<std::valarray<double> > >(const std::valarray<std::valarray<double> >& x)
{
//this implementation uses the OpenMPI MPI_DOUBLE_COMPLEX define
    uint32_t nrof_functions = x.size();
    MPI_Request request;
    if (MPI_Isend(&nrof_functions, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, &request) != 0)//send the number of the transmitted functions to the server
    {
        std::cout<<"An error happened while transmitting the number of functions!"<<std::endl;
        exit(-1);
    }
    std::cout<<"sending to server"<<std::endl;
    double* arr = new double[nrof_functions * x[0].size()];
    //copy into an internal buffer
    for (unsigned int k = 0; k < nrof_functions; ++k)
        for (std::size_t j = 0; j < x[0].size(); ++j)
            arr[k*x[0].size() + j] = x[k][j];
    //send the whole buffer at once
    MPI_Status st;
    MPI_Wait(&request, &st);
    if (MPI_Send(arr, nrof_functions * x[0].size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD) != 0)
    {
        std::cout<<"An error happened while transmitting the actual data!"<<std::endl;
        exit(-1);
    }
    delete [] arr;
}

MPIClient::MPIClient(int argc, char *argv[]) : status_is_uptodate(false)
{
    int numprocs, myid;
    int flag = 0;
    MPI_Initialized(&flag);
    if (!flag) MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    return;
}

MPIClient::~MPIClient()
{
    MPI_Finalize();
}
#endif
