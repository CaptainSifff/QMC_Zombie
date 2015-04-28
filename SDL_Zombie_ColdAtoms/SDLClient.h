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
#include <valarray>
#include "Client.h"
#include "NetError.h"
#include "NetCommonBase.h"
#include "utility.h"
#include "NetClient.h"
#include "charconverters.h"

using namespace std;

class SDLClient : public NetClient
{
public:
    /**
    Constructor call. For the SDLClient we assume that argv[1] is a string that describes the Servers IP
    @param argc argc as given from the main funtion
    @param argv argv as given from the main function
    */
    inline SDLClient(int argc, char *argv[]);
    /**
    With this function you can send something to the Server
    */
    template <typename T>
    inline void sendtoserver(const T&);
    /**
    With this function you can receive data from the server
    */
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
private:
    struct SDL_RAII
    {
        SDL_RAII(){}
        int init(int a)
        {
            return SDL_Init(a);
        }
        ~SDL_RAII()
        {
            SDL_Quit();
        }
    } sdl_raii;
    struct SDLnet_RAII
    {
        SDLnet_RAII(){}
        int init()
        {
            return SDLNet_Init();
        }
        ~SDLnet_RAII()
        {
            SDLNet_Quit();
        }
    } sdlnet_raii;
};

template <typename T>
inline T SDLClient::recvblock()
{
    while (NetClient::isDataReady(1) != 1);
    return recv<T>();
}

template <class T>
struct Priv_SDL_Server
{
    typedef T value_type;
    static inline void sendtoserver(const T& arg, Message& m)
    {
        m.alloc(sizeof(CH<T>));
        CH<T> cd(arg);
        for (unsigned int j = 0; j < cd.size(); ++j)
            m.msg[j] = cd.accessbytes(j);
    }
};

template <class T>
struct Priv_SDL_Server<std::valarray<T> >
{
    typedef std::valarray<T> value_type;
    static inline void sendtoserver(const value_type& arg, Message& m)
    {
        m.alloc(arg.size()*sizeof(CH<T>));
        for (unsigned int k = 0; k < arg.size(); ++k)
        {
            CH<T> cd(arg[k]);
            for (unsigned int j = 0; j < cd.size(); ++j)
                m.msg[k * sizeof(CH<T>) + j] = cd.accessbytes(j);
        }
    }
};

template <class T>
struct Priv_SDL_Server<std::valarray<std::valarray<T> > >
{
    typedef std::valarray<std::valarray<T> > value_type;
    static inline void sendtoserver(const value_type& arg, Message& m)
    {
        const uint32_t nrf = static_cast<uint32_t>(arg.size());
        CH<uint32_t> ui(nrf);
        std::size_t len = arg[0].size();
        m.alloc(sizeof(uint32_t) + nrf * len * sizeof(CH<T>));//the first four bytes hold the number of functions
        for (unsigned int j = 0; j < sizeof(uint32_t); ++j)
            m.msg[j] = ui.accessbytes()[j];
        for (unsigned int i = 0; i < nrf; ++i)
        {
            for (unsigned int k = 0; k < len; ++k)
            {
                CH<T> cd(arg[i][k]);
                for (unsigned int j = 0; j < cd.size(); ++j)
                    m.msg[i *len * sizeof(CH<T>)  + k * sizeof(CH<T>) + j + sizeof(uint32_t)] = cd.accessbytes(j);
            }
        }
    }
};

template <typename T>
inline void SDLClient::sendtoserver(const T& arg)
{
    Message m;
    Priv_SDL_Server<T>::sendtoserver(arg, m);//this prepares the message for sending
    send(m);
    m.free();
    return;
}
//The two Priv_SDL classes are copied over from the server
template <class T>
struct Priv_SDL
{
    typedef T value_type;
    static inline void recvpack(Message& m, value_type& retval)
    {
        CH<T> temp(m.msg);
        retval = temp;
    }
    static inline void sendpack(const T& arg, Message& m)
    {
        m.alloc(sizeof(CH<T>));
        CH<T> cd(arg);
        for (unsigned int j = 0; j < cd.size(); ++j)
            m.msg[j] = cd.accessbytes(j);
    }
};

template <class T>
struct Priv_SDL<std::valarray<T> >
{
    typedef std::valarray<T> value_type;
    static inline void recvpack(Message& m, value_type& retval)
    {
        const std::size_t len = m.len;
        const unsigned int elems = static_cast<unsigned int>(len/sizeof(T));
        retval.resize(elems);
        for (unsigned int k = 0; k < elems; ++k)
        {
            CH<T> temp(&(m.msg[k*sizeof(T)]));
            retval[k] = temp;
        }
    }
    static inline void sendpack(const value_type& arg, Message& m)
    {
        m.alloc(arg.size()*sizeof(CH<T>));
        for (unsigned int k = 0; k < arg.size(); ++k)
        {
            CH<T> cd(arg[k]);
            for (unsigned int j = 0; j < cd.size(); ++j)
                m.msg[k * sizeof(CH<T>) + j] = cd.accessbytes(j);
        }
    }
};


template<typename T>
inline T SDLClient::recv()
{
    Message m = NetClient::recv();
    T ret;
    Priv_SDL<T>::recvpack(m, ret);
    m.free();
    return ret;
}

template<>
inline char SDLClient::recv<char>()
{
    Message m = NetClient::recv();
    char ret = *m.msg;
    m.free();
    return ret;
}

template<>
inline unsigned char SDLClient::recv<unsigned char>()
{
    Message m = NetClient::recv();
    unsigned char ret = *m.msg;
    m.free();
    return ret;
}

/**
This function tries to immediately return
*/
bool SDLClient::isdataready()
{
    return NetClient::isDataReady(0) == 1;
}

SDLClient::SDLClient(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout<<"[SDLClient] No commandline arguments given! Aborting"<<endl;
        exit(-1);
    }
    if (sdl_raii.init(0) == -1)//Initialize SDL
    {
        const string ret("Error in initializing the SDL library!");
        throw SDL_Error(string("[SDLClient] ") + ret);
    }
    if (sdlnet_raii.init() == -1)
    {
        const string ret("Error in initializing the SDL_net library!");
        throw SDL_Error(string("[SDLClient] ") + ret);
    }
    char hostname[256];//single UNIX specification
    gethostname(hostname, 256);//provide sth. like a user name
    uint16_t portBE;
    string hostportspec(argv[1]);
    size_t colonpos = hostportspec.find_first_of(":");
    std::string host;
    if(colonpos == std::string::npos)//port not specified
    {
      SDLNet_Write16(31415, &portBE);
      host = hostportspec;
    }
    else
    {
      uint16_t port;
      host = hostportspec.substr(0, colonpos);
      std::stringstream(hostportspec.substr(colonpos+1)) >> port;
      SDLNet_Write16(port, &portBE);
    }
    uint32_t ip = resolveArgtoIP(host, portBE);
    if (ip == 0xFFFFFFFF)
        throw SDL_Error(string("[SDLClient] Couldn't resolve") + host);
    this->init(ip, portBE, string(hostname));
}
