/***************************************************************************
 *   Copyright (C) 2008,2009 by Florian Goth   *
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
#include <string>
#include "NetCommonBase.h"

class NetClient : protected NetCommonBase<TCPsocket>
{
public:
    inline IPaddress getServerIP();
    inline std::size_t send(const Message&);
    inline Message& recv();
    inline int isDataReady(uint32_t);
    inline NetClient();
    inline ~NetClient();
    inline int init(uint32_t, uint16_t, const string&);
private:
    TCPsocket sock;///< The TCP Socket of the Server
    SDLNet_SocketSet socketset;///< The Socket set where we watch the Server Socket for Activity
    IPaddress serverip;///< The IP of our server
    Message msg;///< a temporaray storage for one received message
    inline Message& getbuf();
    inline void putbuf(Message&);
    string username;
    inline TCPsocket getServerSocket() throw();
};

void NetClient::putbuf(Message& m)
{
    msg = m;
}

Message& NetClient::getbuf()
{
    return msg;
}

NetClient::~NetClient()
{
    SDLNet_FreeSocketSet(socketset);
    SDLNet_TCP_Close(sock);
    return;
}

inline Message& NetClient::recv()
{
    //return NetCommonBase<TCPsocket>::recv(getServerSocket(), data);
    return getbuf();
}

inline std::size_t NetClient::send(const Message& m)
{
    return NetCommonBase<TCPsocket>::send(getServerSocket(), m);
}

/**
@return : 1: Message ready
          0: No Message ready
         -1: An Error occured
*/
int NetClient::isDataReady(uint32_t updateinterval)
{
    uint32_t t1 = SDL_GetTicks();
    int num = SDLNet_CheckSockets(socketset, updateinterval);
    if (num == -1)
    {
        cout<<"SDLNet_CheckSockets: "<<SDLNet_GetError()<<endl;
        return -1;
    }
    if ( (num > 0) && (0 != SDLNet_SocketReady( getServerSocket() )) )
    {
        //There's activity on the Server Port, but we must check if it's a PING msg.
        NetCommonBase<TCPsocket>::recv(getServerSocket(), getbuf() );
        if (getbuf().msg == NULL) return -1;
        if (strncmp(getbuf().msg, "SERVER_PING",  11) == 0)//PING reveived
        {
            static char pong[] = {'C','L', 'I', 'E', 'N','T','_','P','O','N','G', 0};
            NetClient::send(Message(pong, strlen(pong)));
            int dt = SDL_GetTicks() - t1;
            getbuf().free();
            if (dt > static_cast<int>(updateinterval))
                return 0;
            return isDataReady(updateinterval - dt);
        }
        return 1;//it's no Ping
    }
//    cout<<"No Activity"<<endl;
    return 0;
}

inline TCPsocket NetClient::getServerSocket() throw()
{
    return sock;
}

NetClient::NetClient()
{
    return;
}

/**
Initializes the Client
@param ip The IP-Address of the server as a 4-byte integer. Note that the IP-address is Big Endian byte order(aka Network Byte Order).
@param port The port on which to connect as a 2-byte integer. Note that it's in Big-Endian byte order
@return An error code, or a couple of exceptions
*/
int NetClient::init(uint32_t ip, uint16_t port, const string& un)
{
    username = un;
    socketset = SDLNet_AllocSocketSet(1);//Allocate 1 socket for communication with the server
    if (!socketset)
    {
        string ret("SDLNet_AllocSocketSet: " + string(SDLNet_GetError()));
        throw SDL_Error(string("[SDLClient] ") + ret);
    }
    serverip.port = port;
    serverip.host = ip;
    const char* servername = SDLNet_ResolveIP(&serverip);//never free this pointer
    if (servername != NULL)
    {
        cout<<"Connecting to Server: "<<servername<<endl;
    }
    else
        cout<<"Warning: Not able to resolve IP-Address "<<convertIPtoString(ip)<<" !"<<endl;
    unsigned int trys = 0;
    do
    {
        sock = SDLNet_TCP_Open(&serverip);
        if (!sock)
        {
            SDL_Delay(20);//wait a little bit
            ++trys;
        }
    }
    while ((!sock) && (trys < 50));//do multiple tries at getting a connection. In total we wait 50 x 20ms = 1s
    if (!sock)//the server didn't accept the connection
    {
        cout<<"Server was resolved, but didn't accept connection."<<endl;
        cout<<"SDLNet_TCP_Open: "<<SDLNet_GetError()<<endl;
        SDLNet_FreeSocketSet(socketset);
        throw (ErrorConnecting("Server was resolved, but didn't accept connection."));
    }
    if (SDLNet_TCP_AddSocket(socketset, sock) == -1)
    {
        string ret("SDLNet_TCP_AddSocket: " + string(SDLNet_GetError()));
        SDLNet_TCP_Close(sock);
        SDLNet_FreeSocketSet(socketset);
        throw SDL_Error(string("[SDLClient] ") + ret);
    }
    if (send(Message(const_cast<char*>(username.c_str()), username.size() + 1)) == 0)//Send encountered an Error
    {
        string ret("Error while sending User-Name!");
        SDLNet_TCP_Close(sock);
        SDLNet_FreeSocketSet(socketset);
        throw SDL_Error(string("[SDLClient] ") + ret);
    }
    return 0;
}
