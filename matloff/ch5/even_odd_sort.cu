#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

// compare and swap; copies from the f to t, swapping f[i] and
// f[j] if the higher-index value is smaller; it is required that i < j
__device__ void cas(int *f,int *t,int i,int j, int n, int me)
{  
   if (i < 0 || j >= n) return;
   if (me == i) {
      if (f[i] > f[j]) t[me] = f[j];
      else t[me] = f[i];
   } else {  // me == j
      if (f[i] > f[j]) t[me] = f[i];
      else t[me] = f[j];
   }
}

// does one iteration of the sort
__global__ void oekern(int *da, int *daaux, int n, int iter)
{  int bix = blockIdx.x;  // block number within grid
   if (iter % 2) { 
      if (bix % 2) cas(da,daaux,bix-1,bix,n,bix);
      else cas(da,daaux,bix,bix+1,n,bix);
   } else {
      if (bix % 2) cas(da,daaux,bix,bix+1,n,bix);
      else cas(da,daaux,bix-1,bix,n,bix);
   }
}

// sorts the array ha, length n, using odd/even transp. sort; 
// kept simple for illustration, no optimization 
void oddeven(int *ha, int n)
{
   int *da; 
   int dasize = n * sizeof(int);
   cudaMalloc((void **)&da,dasize);
   cudaMemcpy(da,ha,dasize,cudaMemcpyHostToDevice);
   // the array daaux will serve as "scratch space" 
   int *daaux;  
   cudaMalloc((void **)&daaux,dasize);
   dim3 dimGrid(n,1);
   dim3 dimBlock(1,1,1);
   int *tmp;
   for (int iter = 1; iter <= n; iter++) {
      oekern<<<dimGrid,dimBlock>>>(da,daaux,n,iter);
      cudaDeviceSynchronize();
      if (iter < n) {
         // swap pointers
         tmp = da;
         da = daaux;
         daaux = tmp;
      } else
         cudaMemcpy(ha,daaux,dasize,cudaMemcpyDeviceToHost);
   }
}
