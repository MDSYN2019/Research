#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

void m64_init_()
{
 struct rlimit rlim;
 int rv;

 rv=getrlimit(RLIMIT_STACK, &rlim);
 if(rv) {
	perror("getrlimit");
	exit(1);
 }
 rlim.rlim_cur=rlim.rlim_max;
 rv=setrlimit(RLIMIT_STACK, &rlim);
 if(rv) {
	perror("setrlimit");
	exit(1);
 }
}

