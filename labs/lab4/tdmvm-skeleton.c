/* Parallel matrix-vector multiplication with 2D block decomposition
 * The last argument is assumed to be the communicator for a 2D Cartesian
 * topology and we assume the matrix a and vector x are already distribtued
 * with the vector x along the rightmost column of processors.
 */
#include "mpi.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

MatrixVectorMultiply2D(int n, double *a, double *x, double *y, MPI_Comm comm_2d) 
{ 
  int ROW=0, COL=1; /* Improve readability */ 
  int i, j, nlocal; 
  double *py; /* Will store partial dot products */ 
  int npes, dims[2], periods[2], keep_dims[2]; 
  int myrank, mycoords[2]; 
  int source_rank, dest_rank, coords[2]; 
  MPI_Status status; 
  MPI_Comm comm_row, comm_col; 

  /* Get information about the communicator */ 
  MPI_Comm_size(comm_2d, &npes); 
  MPI_Comm_rank(comm_2d, &myrank); 

  /* Compute the size of the square grid */ 
  dims[ROW] = dims[COL] = sqrt(npes); 

  nlocal = n/dims[ROW]; 

  /* Allocate memory for the array that will hold the partial dot-products */ 
  py = malloc(nlocal*sizeof(double)); 

  MPI_Cart_coords(comm_2d, myrank, 2, mycoords); /* Get my coordinates */ 
 
  /* Create the row-based sub-topology */ 
  keep_dims[ROW] = 0; 
  keep_dims[COL] = 1; 
  MPI_Cart_sub(comm_2d, keep_dims, &comm_row); 
 
  /* Create the column-based sub-topology */ 
  /****************************************/


  /****************************************/
 
  /* Redistribute the x vector. */ 
  /* Step 1. The processors along the rightmost column send their data to the diagonal processors */ 
  /* If I'm in the rightmost column but not the last row, send my block
     of the vector to the diagonal processor in my row */ 
  /*****************************************************/





  /*****************************************************/

  /* If I'm on the diagonal but not in the last row, receive the block
     of the vector from the processor in the rightmost column of my row */
  /*****************************************************/





  /*****************************************************/ 
  /* Step 2. Perform a column-wise broadcast with the diagonal process 
             as the root  */ 
  /*******************************************************/


  /*******************************************************/
 
  /* Perform local matrix-vector multiply */ 
  for (i=0; i<nlocal; i++) { 
    py[i] = 0.0; 
    for (j=0; j<nlocal; j++) 
      py[i] += a[i*nlocal+j]*x[j]; 
  } 

  /* Perform the sum-reduction along the rows to add up the partial 
     dot-products and leave the result in the rightmost column */ 
  /********************************************************/


  /********************************************************/
 
  /* free local communicators */
  MPI_Comm_free(&comm_row); /* Free up communicator */ 
  MPI_Comm_free(&comm_col); /* Free up communicator */ 
 
  free(py); 
} 
