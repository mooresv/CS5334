/*
 * exchange.c
 * Uses 2D Cartesian topology to exchange top & bottom rows and left & right columns
 * Should be run with a square number of processes -- e.g., 4, 9, 16, etc.
 */

#include "mpi.h"
#include "math.h"

#define N 8  // Each process has an NxN matrix

int main (int argc, char **argv)
{
	MPI_Status status;
	int nprocs;
	int dims[2];
	int periods[2];
	int reorder = 1;
	MPI_Comm comm2d;
	int my_comm2d_rank;
	int left, right, above, below;
	int matrix[N][N];
	int i, j;
	MPI_Datatype column_type;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

        /* Create the 2D Cartesian communicator comm2d  
	   with wraparound in both dimensions */ 






	// Determine ranks in comm2d of left, right, above, and below neighbors
	


	// Each process initializes its matrix to its comm2d rank
	MPI_Comm_rank(comm2d, &my_comm2d_rank);
	for (i = 0; i < N; i++)
	    for (j = 0; j < N; j++)
		matrix[i][j] = my_comm2d_rank;

	/* Each process exchanges its top and bottom row with its top and bottom
	   neighbor, respectively. */





	/* Each process exchanges its left and right column with its left and right
	 * neighbor, respectively. */
	MPI_Type_vector(N, 1, N, MPI_INT, &column_type);
	MPI_Type_commit(&column_type);




	MPI_Type_free(&column_type);

	MPI_Finalize();
}
