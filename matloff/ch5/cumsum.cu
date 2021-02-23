// for this simple illustration, it is assumed that the code runs in
// just one block, and that the number of threads evenly divides n

// improvements that could be made:
//    1.  change to multiple blocks, to try to use all SMs
//    2.  possibly use shared memory
//    3.  have each thread work on staggered elements of dx, rather than
//        on contiguous ones, to get more efficient bank access

#include <cuda.h>
#include <stdio.h>

__global__ void cumulker(int *dx, int n)
{  
   int me = threadIdx.x;
   int csize  = n / blockDim.x;
   int start = me * csize;
   int i,j,base;
   for (i = 1; i < csize; i++) {
      j = start + i;
      dx[j] = dx[j-1] + dx[j];
   }
   __syncthreads();
   if (me > 0) {
      base = 0;
      for (j = 0; j < me; j++)
         base += dx[(j+1)*csize-1];
   }
   if (me > 0) {
      for (i = start; i < start + csize; i++)
         dx[i] += base;
   }
}


int main(int argc, char **argv)
{
    int n = atoi(argv[1]),  // length of array
        nth = atoi(argv[2]);  // number of threads
    int *ha,  // host array
        *da,  // device array
        nint = n * sizeof(int);
    ha = (int *) malloc(nint);  
    // test example
    for (int i = 0; i < n; i++) ha[i] = i*i % 5;
    if (n < 100) for(int i=0; i<n; i++) printf("%d ",ha[i]);
    printf("\n");
    cudaMalloc((void **)&da,nint);
    cudaMemcpy(da,ha,nint,cudaMemcpyHostToDevice);
    dim3 dimGrid(1,1);  
    dim3 dimBlock(n/nth,1,1);  
    cumulker<<<dimGrid,dimBlock>>>(da,n);
    cudaDeviceSynchronize();
    cudaMemcpy(ha,da,nint,cudaMemcpyDeviceToHost);
    if (n < 100) for(int i=0; i<n; i++) printf("%d ",ha[i]);
    printf("\n");
    free(ha);
    cudaFree(da);
}
