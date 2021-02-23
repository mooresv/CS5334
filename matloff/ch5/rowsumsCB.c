#include <stdio.h>
#include <stdlib.h>
#include <cublas.h>  // required include

int main(int argc, char **argv)
{
   int n = atoi(argv[1]);  // number of matrix rows/cols
   float *hm, // host matrix
         *hrs, // host rowsums vector
         *ones,  // 1s vector for multiply
         *dm, // device matrix
         *drs; // device rowsums vector
   // allocate space on host 
   hm = (float *) malloc(n*n*sizeof(float));
   hrs = (float *) malloc(n*sizeof(float));
   ones = (float *) malloc(n*sizeof(float));
   // as a test, fill hm with consecutive integers, but in column-major
   // order for CUBLAS; also put 1s in ones
   int i,j;
   float t = 0.0;
   for (i = 0; i < n; i++) {
      ones[i] = 1.0;
      for (j = 0; j < n; j++) 
         hm[j*n+i] = t++;
   }
   cublasInit();  // required init
   // set up space on the device
   cublasAlloc(n*n,sizeof(float),(void**)&dm);
   cublasAlloc(n,sizeof(float),(void**)&drs);
   // copy data from host to device
   cublasSetMatrix(n,n,sizeof(float),hm,n,dm,n);
   cublasSetVector(n,sizeof(float),ones,1,drs,1);
   // matrix times vector ("mv")
   cublasSgemv('n',n,n,1.0,dm,n,drs,1,0.0,drs,1);
   // copy result back to host
   cublasGetVector(n,sizeof(float),drs,1,hrs,1);
   // check results
   if (n < 20) for (i = 0; i < n; i++) printf("%f\n",hrs[i]);
   // clean up on device (should call free() on host too)
   cublasFree(dm);
   cublasFree(drs);
   cublasShutdown();
}
