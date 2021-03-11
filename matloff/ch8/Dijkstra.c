// Dijkstra.c

// MPI example program:  Dijkstra shortest-path finder in a
// bidirectional graph; finds the shortest path from vertex 0 to all
// others

// command line arguments:  nv print dbg

// where:  nv is the size of the graph; print is 1 if graph and min
// distances are to be printed out, 0 otherwise; and dbg is 1 or 0, 1
// for debug

// node 0 will both participate in the computation and serve as a
// "manager"

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MYMIN_MSG 0
#define OVRLMIN_MSG 1
#define COLLECT_MSG 2

// global variables (but of course not shared across nodes)

int nv,  // number of vertices
    *notdone, // vertices not checked yet
    nnodes,  // number of MPI nodes in the computation
    chunk,  // number of vertices handled by each node
    startv,endv,  // start, end vertices for this node
    me,  // my node number 
    dbg; 
unsigned largeint,  // max possible unsigned int
         mymin[2],  // mymin[0] is min for my chunk,
                    // mymin[1] is vertex which achieves that min
         othermin[2],  // othermin[0] is min over the other chunks
                       // (used by node 0 only)
                       // othermin[1] is vertex which achieves that min
         overallmin[2],  // overallmin[0] is current min over all nodes,
                         // overallmin[1] is vertex which achieves that min
         *ohd,  // 1-hop distances between vertices; "ohd[i][j]" is
                // ohd[i*nv+j]
         *mind;  // min distances found so far

double T1,T2;  // start and finish times 

void init(int ac, char **av)
{  int i,j,tmp; unsigned u;
   nv = atoi(av[1]);
   dbg = atoi(av[3]);
   MPI_Init(&ac,&av);
   MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
   MPI_Comm_rank(MPI_COMM_WORLD,&me);
   chunk = nv/nnodes;  
   startv = me * chunk; 
   endv = startv + chunk - 1;
   u = -1;
   largeint = u >> 1;
   ohd = malloc(nv*nv*sizeof(int));
   mind = malloc(nv*sizeof(int));
   notdone = malloc(nv*sizeof(int)); 
   // random graph
   // note that this will be generated at all nodes; could generate just
   // at node 0 and then send to others, but faster this way
   srand(9999);
   for (i = 0; i < nv; i++)  
      for (j = i; j < nv; j++)  {
         if (j == i) ohd[i*nv+i] = 0;
         else  {
            ohd[nv*i+j] = rand() % 20;
            ohd[nv*j+i] = ohd[nv*i+j];
         }
      }
   for (i = 0; i < nv; i++)  {
      notdone[i] = 1;
      mind[i] = largeint;
   }
   mind[0] = 0;
   while (dbg) ;  // stalling so can attach debugger
}

// finds closest to 0 among notdone, among startv through endv
void findmymin()
{  int i;
   mymin[0] = largeint; 
   for (i = startv; i <= endv; i++)
      if (notdone[i] && mind[i] < mymin[0])  {
         mymin[0] = mind[i];
         mymin[1] = i;
      }
}

void findoverallmin()
{  int i;
   MPI_Status status;  // describes result of MPI_Recv() call 
   // nodes other than 0 report their mins to node 0, which receives
   // them and updates its value for the global min
   if (me > 0)  
      MPI_Send(mymin,2,MPI_INT,0,MYMIN_MSG,MPI_COMM_WORLD);
   else  {
      // check my own first
      overallmin[0] = mymin[0];
      overallmin[1] = mymin[1];
      // check the others
      for (i = 1; i < nnodes; i++) {
         MPI_Recv(othermin,2,MPI_INT,i,MYMIN_MSG,MPI_COMM_WORLD,&status);
         if (othermin[0] < overallmin[0])  {
            overallmin[0] = othermin[0];
            overallmin[1] = othermin[1];
         }
      }
   }
}

void updatemymind()  // update my mind segment
{  // for each i in [startv,endv], ask whether a shorter path to i
   // exists, through mv
   int i, mv = overallmin[1];
   unsigned md = overallmin[0];
   for (i = startv; i <= endv; i++)
      if (md + ohd[mv*nv+i] < mind[i])
         mind[i] = md + ohd[mv*nv+i];
}

void disseminateoverallmin()
{  int i;
   MPI_Status status;
   if (me == 0)
      for (i = 1; i < nnodes; i++)
         MPI_Send(overallmin,2,MPI_INT,i,OVRLMIN_MSG,MPI_COMM_WORLD);
   else
      MPI_Recv(overallmin,2,MPI_INT,0,OVRLMIN_MSG,MPI_COMM_WORLD,&status);
}

void updateallmind()  // collects all the mind segments at node 0
{  int i;
   MPI_Status status;
   if (me > 0)
      MPI_Send(mind+startv,chunk,MPI_INT,0,COLLECT_MSG,MPI_COMM_WORLD);
   else
      for (i = 1; i < nnodes; i++)
         MPI_Recv(mind+i*chunk,chunk,MPI_INT,i,COLLECT_MSG,MPI_COMM_WORLD,
            &status);
}

void printmind()  // partly for debugging (call from GDB)
{  int i;
   printf("minimum distances:\n");
   for (i = 1; i < nv; i++)
      printf("%u\n",mind[i]);
}

void dowork()
{  int step,  // index for loop of nv steps
       i;
   if (me == 0) T1 = MPI_Wtime();
   for (step = 0; step < nv; step++)  {
      findmymin();
      findoverallmin();
      disseminateoverallmin();
      // mark new vertex as done 
      notdone[overallmin[1]] = 0;  
      updatemymind(startv,endv);
   }
   updateallmind();
   T2 = MPI_Wtime();
}

int main(int ac, char **av)
{  int i,j,print;
   init(ac,av);
   dowork();  
   print = atoi(av[2]);
   if (print && me == 0)  {
      printf("graph weights:\n");
      for (i = 0; i < nv; i++)  {
         for (j = 0; j < nv; j++)  
            printf("%u  ",ohd[nv*i+j]);
         printf("\n");
      }
      printmind();
   }

   if (me == 0) printf("time at node 0: %f\n",(float)(T2-T1));
   MPI_Finalize();
}
