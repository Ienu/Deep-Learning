#include <unistd.h>
#include <iostream>
#include <sys/socket.h>
#include <arpa/inet.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#define MAXSIZE 5

using namespace std;

struct DATA
{
	int x;
	int y;
};

class tcp_server
{
private:
	int socket_fd,accept_fd;
	sockaddr_in myserver;
	sockaddr_in remote_addr;

public:
	tcp_server(int listen_port);
	int recv_msg();
};

tcp_server::tcp_server(int listen_port) 
{
	if(( socket_fd = socket(PF_INET,SOCK_STREAM,IPPROTO_TCP)) < 0 )
	{
		throw "socket() failed";
	}
	memset(&myserver,0,sizeof(myserver));
	myserver.sin_family = AF_INET;
	myserver.sin_addr.s_addr = htonl(INADDR_ANY);
	myserver.sin_port = htons(listen_port);
	if( bind(socket_fd,(sockaddr*) &myserver,sizeof(myserver)) < 0 ) 
	{	
		cout << "ERROR" << endl;
		throw "bind() failed";
	}
	if( listen(socket_fd,10) < 0 ) 
	{
		throw "listen() failed";
	}
}

int tcp_server::recv_msg() 
{
	while(true) 
	{
		socklen_t sin_size = sizeof(struct sockaddr_in);
		if(( accept_fd = accept(socket_fd,(struct sockaddr*) &remote_addr,&sin_size)) == -1 )
		{
			cout << "while error" << endl;
			throw "Accept error!";
			continue;
		}
		printf("Received a connection from %s\n",(char*) inet_ntoa(remote_addr.sin_addr));

		if(!fork())
		{
		  //while(true)
		  //{
				DATA data;
				if( ( read(accept_fd, (char*)&data, sizeof(DATA))) < 0 )
				{
					throw("Read() error!");
					break;
				} 
				else 
				{
				        cout << data.x << " + " << data.y << " = " << data.x + data.y << endl;
				}
				/*if( ( write(accept_fd,buffer,5)) < 0 ) 
				{
				throw("Write() error!");
				}*/
				//}
			close(accept_fd);
		}

	}
	return 0;
}

int main(int argc, char* argv[])
{
	tcp_server ts(atoi(argv[1]));
	ts.recv_msg();
	return 0;
}
