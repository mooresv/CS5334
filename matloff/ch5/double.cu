#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

// CUDA example:  illustrates kernel-allocated shared memory; does
// nothing useful, just copying an array from host to device global,
// then to device shared, doubling it there, then copying back to device
// global then host

__global__ void doubleit(int *dv, int n)
{  extern __shared__ int sv[];
   int me = threadIdx.x;  
   // threads share in copying dv to sv, with each thread copying one
   // element
   sv[me] = dv[me];
   sv[me] = 2 * sv[me];
   dv[me] = sv[me];
}

int main(int argc, char **argv)
{
    int n = atoi(argv[1]);  // number of matrix rows/cols
    int *hv, // host array
        *dv; // device array
    int vsize = n * sizeof(int);  // size of array in bytes
    // allocate space for host array
    hv = (int *) malloc(vsize);  
    // fill test array with consecutive integers
    int t = 0,i;
    for (i = 0; i < n; i++) 
       hv[i] = t++;
    // allocate space for device array 
    cudaMalloc((void **)&dv,vsize);
    // copy host array to device array
    cudaMemcpy(dv,hv,vsize,cudaMemcpyHostToDevice);
    // set up parameters for threads structure
    dim3 dimGrid(1,1);
    dim3 dimBlock(n,1,1);  // all n threads in the same block
    // invoke the kernel; third argument is amount of shared memory
    doubleit<<<dimGrid,dimBlock,vsize>>>(dv,n);
    // wait for kernel to finish
    cudaDeviceSynchronize();
    // copy row array from device to host
    cudaMemcpy(hv,dv,vsize,cudaMemcpyDeviceToHost);
    // check results
    if (n < 10) for(int i=0; i<n; i++) printf("%d\n",hv[i]);
    // clean up
    free(hv);
    cudaFree(dv);
}
