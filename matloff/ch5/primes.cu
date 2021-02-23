#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

// CUDA example:  illustration of shared memory allocation at run time;
// finds primes using classical Sieve of Erathosthenes:  make list of
// numbers 2 to n, then cross out all multiples of 2 (but not 2 itself),
// then all multiples of 3, etc.; whatever is left over is prime; in our
// array, 1 will mean "not crossed out" and 0 will mean "crossed out"

// IMPORTANT NOTE: uses shared memory, in a single block, without
// rotating parts of array in and out of shared memory; thus limited to
// n <= 4000 if have 16K shared memory, and an inefficient use of the
// GPU in any case

// initialize sprimes, 1s for the odds, 0s for the evens; see sieve()
// for the nature of the arguments
__device__ void initsp(int *sprimes, int n, int nth, int me) 
{
   int chunk,startsetsp,endsetsp,val,i;
   sprimes[2] = 1;
   // determine sprimes chunk for this thread to init
   chunk = (n-1) / nth;
   startsetsp = 2 + me*chunk;
   if (me < nth-1) endsetsp = startsetsp + chunk - 1;
   else endsetsp = n;
   // now do the init
   val = startsetsp % 2; 
   for (i = startsetsp; i <= endsetsp; i++) {
      sprimes[i] = val;
      val = 1 - val;
   }
   // make sure sprimes up to date for all
   __syncthreads();
}

// copy sprimes back to device global memory; see sieve() for the nature
// of the arguments
__device__ void cpytoglb(int *dprimes, int *sprimes, int n, int nth, int me) 
{
   int startcpy,endcpy,chunk,i;
   chunk = (n-1) / nth;
   startcpy = 2 + me*chunk;
   if (me < nth-1) endcpy = startcpy + chunk - 1;
   else endcpy = n;
   for (i = startcpy; i <= endcpy; i++) dprimes[i] = sprimes[i];
   __syncthreads();
}

// finds primes from 2 to n, storing the information in dprimes, with
// dprimes[i] being 1 if i is prime, 0 if composite; nth is the number
// of threads (threadDim somehow not recognized)
__global__ void sieve(int *dprimes, int n, int nth)
{
   extern __shared__ int sprimes[];
   int me = threadIdx.x;
   int nth1 = nth - 1;
   // initialize sprimes array, 1s for odds, 0 for evens
   initsp(sprimes,n,nth,me);
   // "cross out" multiples of various numbers m, with each thread doing
   // a chunk of m's; always check first to determine whether m has
   // already been found to be composite; finish when m*m > n
   int maxmult,m,startmult,endmult,chunk,i;
   for (m = 3; m*m <= n; m++) {
      if (sprimes[m] != 0) {
         // find largest multiple of m that is <= n
         maxmult = n / m;
         // now partition 2,3,...,maxmult among the threads
         chunk = (maxmult - 1) / nth;
         startmult = 2 + me*chunk;
         if (me < nth1) endmult = startmult + chunk - 1;
         else endmult = maxmult;
      }
      // OK, cross out my chunk
      for (i = startmult; i <= endmult; i++) sprimes[i*m] = 0;
   }
   __syncthreads();
   // copy back to device global memory for return to host
   cpytoglb(dprimes,sprimes,n,nth,me);
}

int main(int argc, char **argv)
{
    int n = atoi(argv[1]),  // will find primes among 1,...,n
        nth = atoi(argv[2]);  // number of threads
    int *hprimes,  // host primes list
        *dprimes;  // device primes list
    int psize = (n+1) * sizeof(int);  // size of primes lists in bytes
    // allocate space for host list
    hprimes = (int *) malloc(psize);
    // allocate space for device list
    cudaMalloc((void **)&dprimes,psize);
    dim3 dimGrid(1,1);
    dim3 dimBlock(nth,1,1);
    // invoke the kernel, including a request to allocate shared memory
    sieve<<<dimGrid,dimBlock,psize>>>(dprimes,n,nth);
    // check whether we asked for too much shared memory
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) printf("%s\n",cudaGetErrorString(err));
    // wait for kernel to finish
    cudaDeviceSynchronize();
    // copy list from device to host
    cudaMemcpy(hprimes,dprimes,psize,cudaMemcpyDeviceToHost);
    // check results
    if (n <= 1000) for(int i=2; i<=n; i++)
       if (hprimes[i] == 1) printf("%d\n",i);
    // clean up
    free(hprimes);
    cudaFree(dprimes);
}
