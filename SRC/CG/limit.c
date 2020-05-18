#include <stdio.h> // perror
#include <stdlib.h> // exit
#include <sys/time.h> // setrlimit
#include <sys/resource.h> // setrlimit
#include <unistd.h> // setrlimit
void unlimit_stack_(void) {
struct rlimit rlim = { RLIM_INFINITY, RLIM_INFINITY };
if ( setrlimit(RLIMIT_STACK, &rlim) == -1 ) {
perror("setrlimit error");
exit(1);
}
}
