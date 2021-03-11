#include <mpi.h>
#include <stdlib.h>

#define MAX_N 100000
#define MAX_NPROCS 100
#define DATA_MSG 0
#define NEWDATA_MSG 1

int nnodes,  // number of MPI processes
    n,  // size of original array
    me,  // my MPI ID
    has0s[MAX_N],  // original data
    no0s[MAX_N],  // 0-free data
    nno0s;  // number of non-0 elements

int debug;

init(int argc, char **argv)
{  
   int i;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
   MPI_Comm_rank(MPI_COMM_WORLD,&me);
   n = atoi(argv[1]);
   if (me == 0) {
      for (i = 0; i < n; i++) 
         has0s[i] = rand() % 4;
   } else {
      debug = atoi(argv[2]);
      while (debug) ;
   }
}

void managernode()
{  
   MPI_Status status;
   int i;
   int lenchunk;
   lenchunk = n / (nnodes-1);  // assumed divides evenly
   for (i = 1; i < nnodes; i++) {
      MPI_Send(has0s+(i-1)*lenchunk,lenchunk,
         MPI_INT,i,DATA_MSG,MPI_COMM_WORLD);
   }
   int k = 0;
   for (i = 1; i < nnodes; i++) {
      MPI_Recv(no0s+k,MAX_N,
         MPI_INT,i,NEWDATA_MSG,MPI_COMM_WORLD,&status);
      MPI_Get_count(&status,MPI_INT,&lenchunk);
      k += lenchunk;
   }
   nno0s = k;
}

void remov0s(int *oldx, int n, int *newx, int *nnewx)
{  int i,count = 0;
   for (i = 0; i < n; i++)
      if (oldx[i] != 0) newx[count++] = oldx[i];
   *nnewx = count;
}

   if (me == 0 && n < 25) {
      for (i = 0; i < n; i++) printf("%d ",no0s[i]);
      printf("\n");
   }
   MPI_Finalize();
}
