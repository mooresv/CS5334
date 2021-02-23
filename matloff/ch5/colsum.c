#include <stdio.h>
#include <stdlib.h>

// non-CUDA example:  finds col sums of an integer matrix m

// find1elt() finds the colsum of one col of the nxn matrix m, storing the
// result in the corresponding position in the colsum array cs; matrix
// stored as 1-dimensional, row-major order

void find1elt(int *m, int *cs, int n)
{
   int sum=0;
   int topofcol;
   int col,k;
   for (col = 0; col < n; col++) {
      topofcol = col;
      sum = 0;
      for (k = 0; k < n; k++)
         sum += m[topofcol+k*n];
      cs[col] = sum;
   }
}

int main(int argc, char **argv)
{
    int n = atoi(argv[1]);  // number of matrix cols/cols
    int *hm, // host matrix
        *hcs; // host colsums
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
    int cssize = n * sizeof(int);
    hcs = (int *) malloc(cssize);  
    find1elt(hm,hcs,n);
    if (n < 10) for(i=0; i<n; i++) printf("%d\n",hcs[i]);
    // clean up
    free(hm);
    free(hcs);
}
