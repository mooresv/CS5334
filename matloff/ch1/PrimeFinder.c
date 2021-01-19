/* MPI sample program; NOT INTENDED TO BE EFFICIENT as a prime
   finder, either in algorithm or implementation 

   MPI (Message Passing Interface) is a popular package using
   the "message passing" paradigm for communicating between
   processors in parallel applications; as the name implies,
   processors communicate by passing messages using "send" and 
   "receive" functions

   finds and reports the number of primes less than or equal to N

   uses a pipeline approach:  node 0 looks at all the odd numbers (i.e.
   has already done filtering out of multiples of 2) and filters out
   those that are multiples of 3, passing the rest to node 1; node 1
   filters out the multiples of 5, passing the rest to node 2; node 2
   then removes the multiples of 7, and so on; the last node must check
   whatever is left

   note that we should NOT have a node run through all numbers
   before passing them on to the next node, since we would then
   have no parallelism at all; on the other hand, passing on just
   one number at a time isn't efficient either, due to the high
   overhead of sending a message if it is a network (tens of
   microseconds until the first bit reaches the wire, due to
   software delay); thus efficiency would be greatly improved if 
   each node saved up a chunk of numbers before passing them to 
   the next node */

#include <mpi.h>  // mandatory 

#include <stdio.h>
#include <stdlib.h>

#define PIPE_MSG 0  // type of message containing a number to be checked
#define END_MSG 1  // type of message indicating no more data will be coming 

int NNodes,  // number of nodes in computation
    N,  // find all primes from 2 to N 
    Me;  // my node number 
double T1,T2;  // start and finish times 

void Init(int Argc,char **Argv)
{  int DebugWait;  
   N = atoi(Argv[1]); 
   // start debugging section
   DebugWait = atoi(Argv[2]);
   while (DebugWait) ;  // deliberate infinite loop; see below
   /* the above loop is here to synchronize all nodes for debugging;
      if DebugWait is specified as 1 on the mpirun command line, all 
      nodes wait here until the debugging programmer starts GDB at 
      all nodes (via attaching to OS process number), then sets
      some breakpoints, then GDB sets DebugWait to 0 to proceed; */
   // end debugging section
   MPI_Init(&Argc,&Argv);  // mandatory to begin any MPI program 
   // puts the number of nodes in NNodes 
   MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
   // puts the node number of this node in Me 
   MPI_Comm_rank(MPI_COMM_WORLD,&Me); 
   // OK, get started; first record current time in T1 
   if (Me == NNodes-1) T1 = MPI_Wtime();
}

void Node0()
{  int I,ToCheck,Dummy,Error;
   for (I = 1; I <= N/2; I++)  {
      ToCheck = 2 * I + 1;  // latest number to check for div3
      if (ToCheck > N) break;
      if (ToCheck % 3 > 0)  // not divis by 3, so send it down the pipe
         // send the string at ToCheck, consisting of 1 MPI integer, to
         // node 1 among MPI_COMM_WORLD, with a message type PIPE_MSG
         Error = MPI_Send(&ToCheck,1,MPI_INT,1,PIPE_MSG,MPI_COMM_WORLD);
         // error not checked in this code
   }
   // sentinel
   MPI_Send(&Dummy,1,MPI_INT,1,END_MSG,MPI_COMM_WORLD);
}

void NodeBetween()
{  int ToCheck,Dummy,Divisor;
   MPI_Status Status;  
   // first received item gives us our prime divisor
   // receive into Divisor 1 MPI integer from node Me-1, of any message
   // type, and put information about the message in Status
   MPI_Recv(&Divisor,1,MPI_INT,Me-1,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
   while (1)  {
      MPI_Recv(&ToCheck,1,MPI_INT,Me-1,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
      // if the message type was END_MSG, end loop
      if (Status.MPI_TAG == END_MSG) break;
      if (ToCheck % Divisor > 0) 
         MPI_Send(&ToCheck,1,MPI_INT,Me+1,PIPE_MSG,MPI_COMM_WORLD);
   }
   MPI_Send(&Dummy,1,MPI_INT,Me+1,END_MSG,MPI_COMM_WORLD);
}

void NodeEnd()
{  int ToCheck,PrimeCount,I,IsComposite,StartDivisor;
   MPI_Status Status;  
   MPI_Recv(&StartDivisor,1,MPI_INT,Me-1,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
   PrimeCount = Me + 2;  /* must account for the previous primes, which
                            won't be detected below */
   while (1)  {
      MPI_Recv(&ToCheck,1,MPI_INT,Me-1,MPI_ANY_TAG,MPI_COMM_WORLD,&Status);
      if (Status.MPI_TAG == END_MSG) break;
      IsComposite = 0;
      for (I = StartDivisor; I*I <= ToCheck; I += 2)
         if (ToCheck % I == 0)  {
            IsComposite = 1;
            break;
         }
      if (!IsComposite) PrimeCount++;  
   }
   /* check the time again, and subtract to find run time */
   T2 = MPI_Wtime();
   printf("elapsed time = %f\n",(float)(T2-T1));
   /* print results */
   printf("number of primes = %d\n",PrimeCount);
}

int main(int argc,char **argv)
{  Init(argc,argv);
   // all nodes run this same program, but different nodes take
   // different actions
   if (Me == 0) Node0();
   else if (Me == NNodes-1) NodeEnd();
        else NodeBetween();
   // mandatory for all MPI programs 
   MPI_Finalize();
}

/* explanation of "number of items" and "status" arguments at the end 
   of MPI_Recv():

   when receiving a message you must anticipate the longest possible
   message, but the actual received message may be much shorter than
   this; you can call the MPI_Get_count() function on the status
   argument to find out how many items were actually received

   the status argument will be a pointer to a struct, containing the
   node number, message type and error status of the received
   message

   say our last parameter is Status; then Status.MPI_SOURCE
   will contain the number of the sending node, and 
   Status.MPI_TAG will contain the message type; these are
   important if used MPI_ANY_SOURCE or MPI_ANY_TAG in our
   node or tag fields but still have to know who sent the
   message or what kind it is */
