#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

// CUDA example:  finds row sums of an integer matrix m

// find1elt() finds the rowsum of one row of the nxn matrix m, storing the
// result in the corresponding position in the rowsum array rs; matrix
// stored as 1-dimensional, row-major order

__global__ void find1elt(int *m, int *rs, int n)
{
   int rownum = blockIdx.x;  // this thread will handle row # rownum
   int sum = 0;
   for (int k = 0; k < n; k++)
      sum += m[rownum*n+k];
   rs[rownum] = sum;
}

int main(int argc, char **argv)
{
    int n = atoi(argv[1]);  // number of matrix rows/cols
    int *hm, // host matrix
        *dm, // device matrix
        *hrs, // host rowsums
        *drs; // device rowsums
    int msize = n * n * sizeof(int);  // size of matrix in bytes
    // allocate space for host matrix
    hm = (int *) malloc(msize);  
    // as a test, fill matrix with consecutive integers
    int t = 0,i,j;
    for (i = 0; i < n; i++) {
       for (j = 0; j < n; j++) {
          hm[i*n+j] = t++;
       }
    }
    // allocate space for device matrix 
    cudaMalloc((void **)&dm,msize);
    // copy host matrix to device matrix
    cudaMemcpy(dm,hm,msize,cudaMemcpyHostToDevice);
    // allocate host, device rowsum arrays
    int rssize = n * sizeof(int);
    hrs = (int *) malloc(rssize);  
    cudaMalloc((void **)&drs,rssize);
    // set up parameters for threads structure
    dim3 dimGrid(n,1);  // n blocks 
    dim3 dimBlock(1,1,1);  // 1 thread per block
    // invoke the kernel
    find1elt<<<dimGrid,dimBlock>>>(dm,drs,n);
    // wait for kernel to finish
    cudaDeviceSynchronize();
    // copy row vector from device to host
    cudaMemcpy(hrs,drs,rssize,cudaMemcpyDeviceToHost);
    // check results
    if (n < 10) for(int i=0; i<n; i++) printf("%d\n",hrs[i]);
    // clean up
    free(hm);
    cudaFree(dm);
    free(hrs);
    cudaFree(drs);
}
