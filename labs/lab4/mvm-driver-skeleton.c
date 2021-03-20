/* test driver program for 1D matrix-vector multiplication */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

extern int MatrixVectorMultiply(int n, double *a, double *b, double *c, MPI_Comm comm);

int main (int argc, char *argv[])
{
int   numtasks, taskid;
int i, j, N, n, nlocal; 
double *a, *x, *y, *ycheck, *alocal, *xlocal, *ylocal;
int proc;
MPI_Status status;
MPI_Request request;

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

if (argc != 2) {
   if (taskid == 0) {
     fprintf(stderr, "Usage: %s <n>\n", argv[0]);
     fprintf(stderr, "where n is a multiple of the number of tasks.\n");
   }   
   MPI_Finalize();
   exit(0);
}

/* Read row/column dimenstion from command line */
n = atoi(argv[1]);

N=n*n; //size of matrix n x n
nlocal = n/numtasks;

if (taskid == 0) {
  a = (double *) malloc(N*sizeof(double));
  x = (double *) malloc(n*sizeof(double));
  y = (double *) malloc(n*sizeof(double));
  ycheck = (double *) malloc(n*sizeof(double));

  /* Initialize a and x */
  for (i = 0;i < n;i++) {
    for (j = 0;j < n; j++)
      a[i*n +j] = 2*i+j;
    x[i] = i;
  }
}

alocal = (double *) malloc(n*nlocal*sizeof(double));
xlocal = (double *) malloc(nlocal*sizeof(double));
ylocal = (double *) malloc(nlocal*sizeof(double));

/* Distribute a and x in 1D row distribution */
/*********************************************/


/*********************************************/

/* Each process calls MatrixVectorMultiply */

MatrixVectorMultiply(n, alocal, xlocal, ylocal, MPI_COMM_WORLD);

/* Gather results back to root process */
/***************************************/

/***************************************/

/* Check results */

if (taskid == 0) { 
  for (i = 0; i<n; i++) {
     ycheck[i] = 0;
     for (j=0; j<n; j++) 
	ycheck[i] += a[i*n+j]*x[j];
     if (ycheck[i] != y[i])
        printf("discrepancy: ycheck[%d]=%f, y[%d]=%f\n", i, ycheck[i], i, y[i]);
  }
  printf("Done with mvm, y[%d] = %f\n", n-1, y[n-1]);
}

MPI_Finalize();

}
