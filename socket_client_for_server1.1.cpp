#include <stdio.h>
#include <Winsock2.h>

#pragma comment(lib,"WS2_32.lib") 

struct DATA
{
	int x;
	int y;
};

void main()
{
	WORD wVersionRequested;
	WSADATA wsaData;
	int err;
	
	wVersionRequested = MAKEWORD( 1, 1 );
	
	err = WSAStartup( wVersionRequested, &wsaData );
	if ( err != 0 ) {
		return;
	}
	
	if ( LOBYTE( wsaData.wVersion ) != 1 ||
        HIBYTE( wsaData.wVersion ) != 1 ) {
		WSACleanup( );
		return;
	}
	SOCKET sockClient=socket(AF_INET,SOCK_STREAM,0);
	
	SOCKADDR_IN addrSrv;
	addrSrv.sin_addr.S_un.S_addr=inet_addr("192.168.1.23");
	addrSrv.sin_family=AF_INET;
	addrSrv.sin_port=htons(30000);
	connect(sockClient,(SOCKADDR*)&addrSrv,sizeof(SOCKADDR));

	DATA data;
	data.x = 2;
	data.y = 2;

	send(sockClient, (char*)(&data), sizeof(DATA), 0);
	//char recvBuf[50];
	//recv(sockClient,recvBuf,50,0);
	//printf("%s\n",recvBuf);
	
	closesocket(sockClient);
	WSACleanup();
}
