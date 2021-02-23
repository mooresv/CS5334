#include<cuda.h>
#include <stdio.h>

int main()
{
   cudaDeviceProp Props;
   cudaGetDeviceProperties( &Props,0);

   printf("shared mem: %d)\n", Props.sharedMemPerBlock);
   printf("max threads/block: %d\n",Props.maxThreadsPerBlock);
   printf("max blocks: %d\n",Props.maxGridSize[0]);
   printf("total Const mem: %d\n",Props.totalConstMem);
}
