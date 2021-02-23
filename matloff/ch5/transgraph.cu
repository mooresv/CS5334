// takes a graph adjacency matrix for a directed graph, and converts it
// to a 2-column matrix of pairs (i,j), meaning an edge from vertex i to
// vertex j; the output matrix must be in lexicographical order

// not claimed efficient, either in speed or in memory usage

#include <cuda.h>
#include <stdio.h>

// needs -lrt link flag for C++
#include <time.h>
float timediff(struct timespec t1, struct timespec t2)
{  if (t1.tv_nsec > t2.tv_nsec) {
         t2.tv_sec -= 1;
         t2.tv_nsec += 1000000000;
   }
   return t2.tv_sec-t1.tv_sec + 0.000000001 * (t2.tv_nsec-t1.tv_nsec);
}

// tgkernel1() finds the number of 1s to be handled by a thread, used
// to determine where in the output matrix a thread writes its portion

// arguments:
//    dadjm:  the adjacency matrix (NOT assumed symmetric), 1 for edge, 0
//            otherwise; note: matrix is overwritten by the function
//    n:  number of rows and columns of adjm
//    dcounts:  output array, counts of 1s

__global__ void tgkernel1(int *dadjm, int n, int *dcounts)
{  int tot1s,j;
   // need to find my thread number among the totality of all
   // threads in all blocks
   int me = blockDim.x * blockIdx.x + threadIdx.x;
   tot1s = 0;
   for (j = 0; j < n; j++) {
      if (dadjm[n*me+j] == 1) {
         dadjm[n*me+tot1s++] = j;
      }
   dcounts[me] = tot1s;
   }
}

// tgkernel2() has the given thread write its rows into the output
// matrix

__global__ void tgkernel2(int *dadjm, int n, 
   int *dcounts, int *dstarts, int *doutm)
{  int outrow,num1si,j;
   int me = blockDim.x * blockIdx.x + threadIdx.x;
   // fill in this thread's portion of doutm
   outrow = dstarts[me];
   num1si = dcounts[me];
   if (num1si > 0) {
      for (j = 0; j < num1si; j++) {
         doutm[2*outrow+2*j] = me;
         doutm[2*outrow+2*j+1] = dadjm[n*me+j];
     }
   }
}

// replaces counts by cumulative counts
void cumulcounts(int *c, int *s, int n)
{  int i;
   s[0] = 0;
   for (i = 1; i < n; i++) {
      s[i] = s[i-1] + c[i-1];
   }
}

int *transgraph(int *hadjm, int n, int *nout, int gsize, int bsize)
{  int *dadjm;  // device adjacency matrix
   int *houtm;  // host output matrix
   int *doutm;  // device output matrix
   int *hcounts;  // host counts vector
   int *dcounts;  // device counts vector
   int *hstarts;  // host starts vector
   int *dstarts;  // device starts vector
   hcounts = (int *) malloc(n*sizeof(int));
   hstarts = (int *) malloc(n*sizeof(int));
   cudaMalloc((void **)&dadjm,n*n*sizeof(int));
   cudaMalloc((void **)&dcounts,n*sizeof(int));
   cudaMalloc((void **)&dstarts,n*sizeof(int));
   houtm = (int *) malloc(n*n*sizeof(int));
   cudaMalloc((void **)&doutm,n*n*sizeof(int));
   cudaMemcpy(dadjm,hadjm,n*n*sizeof(int),cudaMemcpyHostToDevice);
   dim3 dimGrid(gsize,1);
   dim3 dimBlock(bsize,1,1);
   // calculate counts and starts first
   tgkernel1<<<dimGrid,dimBlock>>>(dadjm,n,dcounts);
   cudaMemcpy(hcounts,dcounts,n*sizeof(int),cudaMemcpyDeviceToHost);
   cumulcounts(hcounts,hstarts,n);
   *nout = hstarts[n-1] + hcounts[n-1];
   cudaMemcpy(dstarts,hstarts,n*sizeof(int),cudaMemcpyHostToDevice);
   tgkernel2<<<dimGrid,dimBlock>>>(dadjm,n,dcounts,dstarts,doutm);
   cudaMemcpy(houtm,doutm,2*(*nout)*sizeof(int),cudaMemcpyDeviceToHost);
   free(hcounts);
   free(hstarts);
   cudaFree(dadjm);
   cudaFree(dcounts);
   cudaFree(dstarts);
   return houtm;
}

int main(int argc, char **argv)
{  int i,j;
   int *adjm;  // host adjacency matrix
   int *outm;  // host output matrix
   int n = atoi(argv[1]);
   int gsize = atoi(argv[2]);
   int bsize = atoi(argv[3]);
   int nout;
   adjm = (int *) malloc(n*n*sizeof(int));
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         if (i == j) adjm[n*i+j] = 0;
         else adjm[n*i+j] = rand() % 2;
   if (n < 10) {
      printf("adjacency matrix: \n");
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) printf("%d ",adjm[n*i+j]);
         printf("\n");
      }
   }

   struct timespec bgn,nd;
   clock_gettime(CLOCK_REALTIME, &bgn);

   outm = transgraph(adjm,n,&nout,gsize,bsize);
   printf("num rows in out matrix = %d\n",nout);
   if (nout < 50) {
      printf("out matrix: \n");
      for (i = 0; i < nout; i++)
         printf("%d %d\n",outm[2*i],outm[2*i+1]);
   }

   clock_gettime(CLOCK_REALTIME, &nd);
   printf("%f\n",timediff(bgn,nd));
}
