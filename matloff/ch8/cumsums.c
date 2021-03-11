// finds cumulative sums in the array x

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_N 10000000  
#define MAX_NODES 10

int nnodes,  // number of MPI processes
    n,  // size of x
    me,  // MPI rank of this node
    // full data for node 0, part for the rest
    x[MAX_N],  
    csums[MAX_N],  // cumulative sums for this node
    maxvals[MAX_NODES];  // the max values at the various nodes 

int debug; 

void init(int argc, char **argv)
{  
   int i;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
   MPI_Comm_rank(MPI_COMM_WORLD,&me); 
   n = atoi(argv[1]); 
   // test data
   if (me == 0) {
      for (i = 0; i < n; i++) 
         x[i] = rand() % 32;
   } 
   debug = atoi(argv[2]); 
   while (debug) ;
}

void cumulsums()
{  
   MPI_Status status;
   int i,lenchunk,sum,node; 
   lenchunk = n / nnodes;  // assumed to divide evenly
   // note that node 0 will participate in the computation too
   MPI_Scatter(x,lenchunk,MPI_INT,x,lenchunk,MPI_INT,
      0,MPI_COMM_WORLD);
   sum = 0;
   for (i = 0; i < lenchunk; i++) {
      csums[i] = sum + x[i];
      sum += x[i];
   }
   MPI_Gather(&csums[lenchunk-1],1,MPI_INT,
      maxvals,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(maxvals,nnodes,MPI_INT,0,MPI_COMM_WORLD);
   if (me > 0) {
      sum = 0;
      for (node = 0; node < me; node++) {
         sum += maxvals[node];
      }
      for (i = 0; i < lenchunk; i++) 
         csums[i] += sum;
   }
   MPI_Gather(csums,lenchunk,MPI_INT,csums,lenchunk,MPI_INT,
      0,MPI_COMM_WORLD);
}

int main(int argc,char **argv)
{  
   int i;
   init(argc,argv);
   if (me == 0 && n < 25) {
      for (i = 0; i < n; i++) printf("%d ",x[i]);
      printf("\n");
   }
   cumulsums();
   if (me == 0 && n < 25) {
      for (i = 0; i < n; i++) printf("%d ",csums[i]);
      printf("\n");
   }
   MPI_Finalize();
}
