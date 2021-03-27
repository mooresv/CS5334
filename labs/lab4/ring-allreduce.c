//
// perform all-reduce of data stored on each process by rotating each piece of data all the
// way round the ring. At each iteration, a process receives some data from the left, adds the 
// value to its running total, then passes the data it has just received on to the right.
//
// Use virtual topologies, here cartesian 1D
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define TRUE  1
#define FALSE 0
#define NDIMS 1    

int main(){
  MPI_Comm comm1d;
  int direction,disp;
  int dims[NDIMS];
  int period[NDIMS];
  int coords[NDIMS];
  int reorder;
  int comm_world_rank, comm1d_rank, comm1d_cart_rank, comm_size, tag;
  int i;
  int left, right, addon, passon, sum;
  MPI_Status status;
  MPI_Request request;

  tag =1;

  MPI_Init(NULL,NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size); 

  // Cartesian topology
  dims[0] = 0;
  period[0] = TRUE;    // wraparound ring
  reorder = TRUE;     
  direction = 0;       // shift along the first index
  disp = 1;            // shift by 1

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_world_rank);
  MPI_Dims_create(comm_size,NDIMS,dims);
  MPI_Cart_create(MPI_COMM_WORLD,NDIMS,dims,period,reorder,&comm1d);
  MPI_Comm_rank(comm1d,&comm1d_rank);
  MPI_Cart_coords(comm1d, comm1d_rank, 1, coords);
  MPI_Cart_rank(comm1d, coords, &comm1d_cart_rank);
  MPI_Cart_shift(comm1d,direction,disp,&left,&right);

  sum = 0;

  // Initialise local values to:
  passon = (coords[0]+1)*(coords[0]+1);

  // Use non-blocking point-to-point communication

  for(i=0;i<comm_size;i++){
    MPI_Issend(&passon,1,MPI_INT,right,tag,comm1d,&request);
    MPI_Recv(&addon,1,MPI_INT,left,tag,comm1d,&status);
    MPI_Wait(&request,&status);

    sum = sum + addon;
    passon = addon;

  }


  printf( "The sum is: %d on comm_world_rank %d, comm1d_rank %d comm1d_cart_rank %d, comm1d_coord %d\n",sum, comm_world_rank, comm1d_rank, comm1d_cart_rank, coords[0]);
  MPI_Finalize();

}
