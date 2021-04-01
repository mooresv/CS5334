/**************************************
      Conway Game of Life
 2-process domain decomposition;
 domain decomposed with horizontal
 line, i.e., top half and bottom half
***************************************/
 
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NI 200
#define NJ 200
#define NSTEPS 500
#define REORDER 1
void main(int argc, char *argv[]){
  int i, j, n, im, ip, jm, jp, nsum, isum, myisum, isumloc, nprocs ,myid, my2drank;
  int N;
  int ig, jg, i1g, i2g, j1g, j2g, ninom, njnom, ninj, i1, i2, i2m, j1, j2, j2m, ni, nj;
  int niproc, njproc;  /* no. procs in each direction */
  int **old, **new, *old1d, *new1d;
  int dims[2], mycoords[2], periods[2], coords[2];
  int left, right, above, below, upper_left, upper_right, lower_left, lower_right;
  MPI_Comm comm_2d;
  MPI_Datatype column_type;
  MPI_Status status;
  float x;
 
  /* initialize MPI */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  /* nominal number of points per proc. in each direction
    (without ghost cells; assume numbers divide evenly) */ 
  niproc = sqrt(nprocs);  /* divide domain in i and j directions */
  njproc = niproc;
  ninom = NI/niproc;
  njnom = NJ/njproc;
 
  dims[0] = niproc;
  dims[1] = njproc;
  periods[0] = 1;
  periods[1] = 1;

  /* create 2d Cartesian communicator, find coords and neighbor ranks */













  //printf("my2drank: %d, left: %d, right: %d, above: %d, below: %d, upper_left: %d, upper_right: %d, lower_left: %d, lower_right: %d\n", my2drank, left, right, above, below, upper_left, upper_right, lower_left, lower_right);

  /* domain decomposition */
 
  /* global starting and ending indices (without ghost cells) */
  i1g = (mycoords[0]*ninom) + 1;
  i2g = i1g + ninom - 1;
  j1g = (mycoords[1]*njnom) + 1;
  j2g = j1g + ninom - 1;
 
  /* local starting and ending indices, including ghost cells */
  i1  = 0;
  i2  = ninom + 1;
  i2m = i2-1;
  j1  = 0;
  j2  = njnom+1;
  j2m = j2-1;
 
  /* allocate arrays; want elements to be contiguous,
     so allocate 1-D arrays, then set pointer to each row
     (old and new) to allow use of array notation for convenience */
  ni = i2-i1+1;
  nj = j2-j1+1;
  ninj = ni*nj;
  old1d = malloc(ninj*sizeof(int));
  new1d = malloc(ninj*sizeof(int));
  old   = malloc(ni*sizeof(int*));
  new   = malloc(ni*sizeof(int*));
  for(i=0; i<ni; i++){
    old[i] = &old1d[i*nj];
    new[i] = &new1d[i*nj];
  }

  /*  Initialize elements of old to 0 or 1.
      We're doing some sleight of hand here to make sure we
      initialize to the same values as in the serial code.
      The rand() function is called for every i and j, even
      if they are not on the current process, to get the same
      random distribution as the serial case, but they are
      only used if this (i,j) resides on the current process. */ 
 
   srand(100000);
   for(ig=1; ig<=NI; ig++){
    for(jg=1; jg<=NJ; jg++){
      x = rand()/((float)RAND_MAX + 1);
 
      // if this i is on the current process
      if( ig >= i1g && ig <= i2g && jg >= j1g && jg <= j2g ){
        // local i and j indices, accounting for lower ghost cell 
        i = ig - i1g + 1;
        j = jg - j1g + 1;
        if(x<0.5){
	  old[i][j] = 0;
	} else {
          old[i][j] = 1;
	}
     }
    }
  }
  //if (my2drank == 0) printf("Initial matrix\n");
  myisum = 0;
  for(i=1; i<i2; i++){
    for(j=1; j<j2; j++){
      myisum = myisum + old[i][j];
   //   if (my2drank == 0) printf("%d ", old[i][j]);
    }
    //  if (my2drank == 0) printf("\n");
  }
  /* Print initial number of live cells. */ 
  MPI_Reduce(&myisum, &isum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(myid == 0) printf("Initial number of live cells = %d\n", isum);

  /* Create a derived type for a column of a processor's block.
   * There are ninom blocks, each containing 1 element, with a stride of
   * nj between blocks */


  /*  Iterate */
  for(n=0; n<NSTEPS; n++){
   /* transfer data to ghost cells */
    if(nprocs == 1){
      for(i=1; i<i2; i++){          /* left and right columns */
        old[i][0]  = old[i][j2m];
        old[i][j2] = old[i][1];
      }
      for(j=1; j<j2; j++){          /* top and bottom rows */
        old[0][j]  = old[i2m][j];
        old[i2][j] = old[1][j];
      }
      old[0][0]   = old[i2m][j2m];  /* corners */
      old[0][j2]  = old[i2m][1];
      old[i2][0]  = old[1][j2m];
      old[i2][j2] = old[1][1];
   } else{
      /* top and bottom rows */


      /* left and right columns using column_type derived datatype */


      /* corners */



   }
    for(i=1; i<i2; i++){
      for(j=1; j<j2; j++){     
        im = i-1;
        ip = i+1;
        jm = j-1;
        jp = j+1;
        nsum =  old[im][jp] + old[i][jp] + old[ip][jp]
                  + old[im][j ]              + old[ip][j ] 
                  + old[im][jm] + old[i][jm] + old[ip][jm];
        switch(nsum){
        case 3:
          new[i][j] = 1;
          break;
        case 2:
          new[i][j] = old[i][j];
          break;
        default:
          new[i][j] = 0;
        }
      }
    }
 
    /* copy new state into old state */    
    for(i=1; i<i2; i++){
      for(j=1; j<j2; j++){
        old[i][j] = new[i][j];
      }
    }
  }
 
  /*  Iterations are done; sum the number of live cells */
  myisum = 0;
  for(i=1; i<i2; i++){
    for(j=1; j<j2; j++){
      myisum = myisum + new[i][j];
    }
  }
 
  /* Print final number of live cells. */ 
  MPI_Reduce(&myisum, &isum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(myid == 0) printf("Final number of live cells = %d\n", isum);
  MPI_Finalize();
}
