#include<stdlib.h>
#include<stdio.h>
 
#define N 2048
typedef long long long_int;
 
int main()
{
 
  int dims[3] = {N,N,N};
 
  long_int size;
  double *array;
  int i;
  long_int li;
 
  size = 1;
  for (i=0; i < 3; i++)
    size *= dims[i];
  printf("Size %lld\n", size);
 
  array = new double[size];
 
  for (li=0; li < size; li++)
    array[li] = (double) (li);
 
  li = size - 1;
  printf("Last %16.12e\n", array[li]);
}
