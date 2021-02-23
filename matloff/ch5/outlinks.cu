#include <cuda.h>
#include <stdio.h>

// CUDA example:  finds mean number of mutual outlinks, among all pairs
// of Web sites in our set; in checking all (i,j) pairs, thread k will
// handle all i such that i mod totth = k, where totth is the number of
// threads

// procpairs() processes all pairs for a given thread
__global__ void procpairs(int *m, int *tot, int n)
{  int totth = gridDim.x * blockDim.x,  // total number of threads
       // need to find my thread number among the totality of all
       // threads in all blocks
       me = blockIdx.x * blockDim.x + threadIdx.x;  
   int i,j,k,sum = 0; 
   for (i = me; i < n; i += totth) {  // do various rows i
      for (j = i+1; j < n; j++) {  // do all rows j > i
         for (k = 0; k < n; k++)
            sum += m[n*i+k] * m[n*j+k];
      }
   }
   atomicAdd(tot,sum);
}

int main(int argc, char **argv)
{  int n = atoi(argv[1]),  // number of vertices
       nblk = atoi(argv[2]);  // number of blocks
    int *hm, // host matrix
        *dm, // device matrix
        htot, // host grand total
        *dtot; // device grand total
    int msize = n * n * sizeof(int);  // size of matrix in bytes
    // allocate space for host matrix
    hm = (int *) malloc(msize);  
    // as a test, fill matrix with random 1s and 0s
    int i,j;
    for (i = 0; i < n; i++) {
       hm[n*i+i] = 0;
       for (j = 0; j < n; j++) {
          if (j != i) hm[i*n+j] = rand() % 2;
       }
    }
    // allocate space for device matrix 
    cudaMalloc((void **)&dm,msize);
    // copy host matrix to device matrix
    cudaMemcpy(dm,hm,msize,cudaMemcpyHostToDevice);
    htot = 0;
    // set up device total and initialize it
    cudaMalloc((void **)&dtot,sizeof(int));
    cudaMemcpy(dtot,&htot,sizeof(int),cudaMemcpyHostToDevice);
    // set up parameters for threads structure
    dim3 dimGrid(nblk,1);
    dim3 dimBlock(192,1,1);
    // invoke the kernel
    procpairs<<<dimGrid,dimBlock>>>(dm,dtot,n);
    // wait for kernel to finish
    cudaDeviceSynchronize();
    // copy total from device to host
    cudaMemcpy(&htot,dtot,sizeof(int),cudaMemcpyDeviceToHost);
    // check results
    if (n <= 15) {
       for (i = 0; i < n; i++) {
          for (j = 0; j < n; j++) 
             printf("%d ",hm[n*i+j]);
          printf("\n");
       }
    }
    printf("mean = %f\n",htot/float((n*(n-1))/2));
    // clean up
    free(hm);
    cudaFree(dm);
    cudaFree(dtot);
}
