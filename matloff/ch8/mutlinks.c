// MPI solution to the mutual outlinks problem

// adjacency matrix m is global at each node, broadcast from node 0

// assumes m is nxn, and number of nodes is < n

// for each node i, check all possible pairing nodes j > i; the various
// nodes work on values of i in a Round Robin fashion, with node k
// handling all i for which i mod nnodes = k

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#define MAXLENGTH 10000000  

int nnodes,  // number of MPI processes
    n,  // size of x
    me,  // MPI rank of this node
    m[MAXLENGTH],  // adjacency matrix
    grandtot;  // grand total of all counts of mutuality

// get adjacency matrix, in this case just by simulation
void getm() 
{  int i;
   for (i = 0; i < n*n; i++) 
      m[i] = rand() % 2;
}

void init(int argc, char **argv)
{  
   int i;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
   MPI_Comm_rank(MPI_COMM_WORLD,&me); 
   n = atoi(argv[1]); 
   if (me == 0) {
      getm();  // get the data (app-specific)
   } 
}

// convert 2-D subscript to 1-D
int twod2oned(int n,int i,int j) 
{  return n * i + j;  }

void mutlinks()
{  
   int i,j,k,tot; 
   MPI_Bcast(m,n*n,MPI_INT,0,MPI_COMM_WORLD);
   tot = 0;
   for (i = me; i < n-1; i += nnodes) {
      for (j = i+1; j < n; j++) {
         for (k = 0; k < n; k++) 
            tot += m[twod2oned(n,i,k)] * m[twod2oned(n,j,k)];
      }
   }
   MPI_Reduce(&tot,&grandtot,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
}

int main(int argc,char **argv)
{  int i,j;
   init(argc,argv);
   if (me == 0 && n < 5) {  // check test input
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) printf("%d ",m[twod2oned(n,i,j)]);
         printf("\n");
      }
   }
   mutlinks();
   if (me == 0) printf("%f\n",((float) grandtot)/(n*(n-1)/2));
   MPI_Finalize();
}
