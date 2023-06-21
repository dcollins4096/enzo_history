#include<stdio.h>
#include<sys/time.h>
#include<sys/resource.h>
#include<unistd.h>
#include<assert.h>

long long int mused(void)
{

#ifdef MEM_TRACE
  struct rusage temp;
  long long int bytes;
  int result;

  result = getrusage(RUSAGE_SELF, &temp);
  if( result == 0 ) {
    bytes = ((long long int) (1024)) * ((long long int) temp.ru_maxrss);
  } else {
    bytes = ((long long int) (0));
  }
  return(bytes);
#else
  return(0);
#endif

}
