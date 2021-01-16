// PrimesThreads.c

// threads-based program to find the number of primes between 2 and n;
// uses the Sieve of Eratosthenes, deleting all multiples of 2, all
// multiples of 3, all multiples of 5, etc. 

// for illustration purposes only; NOT claimed to be efficient

// usage:  primesthreads n num_threads

#include "PrimeThreads.h"

// "crosses out" all odd multiples of k
void crossout(int k)
{  int i;
   for (i = 3; i*k <= n; i += 2)  {
      prime[i*k] = 0;
   }
}

// each thread runs this routine
void *worker(void *tn)  // tn is the thread number (0,1,...)
{  int lim,base;
   long work = 0;  // amount of work done by this thread
   // no need to check multipliers bigger than sqrt(n)
   lim = sqrt(n);
   do  {
      // get next sieve multiplier, avoiding duplication across threads
      // lock the lock
      pthread_mutex_lock(&nextbaselock);
      base = nextbase;
      nextbase += 2;
      // unlock
      pthread_mutex_unlock(&nextbaselock);
      if (base <= lim)  {
         // don't bother crossing out if base known composite
         if (prime[base])  {
            crossout(base);
            work++;  // log work done by this thread
         }
      }
      else pthread_exit((void*) work);
      //else return work; 
   } while (1);
}

int main(int argc, char **argv)
{  int nprimes;  // number of primes found 
   void *status;
   long t;
   int i;
   pthread_attr_t attr;
   pthread_t id[MAX_THREADS];
   long totalwork = 0;

   n = atoi(argv[1]);
   nthreads = atoi(argv[2]);
   // mark all even numbers nonprime, and the rest "prime until
   // shown otherwise"
   for (i = 3; i <= n; i++)  {
      if (i%2 == 0) prime[i] = 0;
      else prime[i] = 1;
   }
   nextbase = 3;

   /* Initialize and set thread detached attribute */
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

   // get threads started
   for (t = 0; t < nthreads; t++)  {
      // this call says create a thread, record its ID in the array
      // id, and get the thread started executing the function worker(), 
      // passing the argument i to that function
      //printf("calling pthread_create, t=%ld\n", t);
      pthread_create(&id[t],&attr,worker, (void *) t);
   }

   // wait for all done
   //printf("nthreads: %d\n", nthreads);
   for (t = 0; t < nthreads; t++)  {
      // this call says wait until thread number id[i] finishes
      // execution, and to assign the return value of that thread to our
      // local variable work here
      //printf("calling pthread_join, t=%ld\n", t);
      pthread_join(id[t],&status);
      totalwork += (long) status;
      printf("%ld values of base done\n",(long)status);
   }

   // report results
   nprimes = 1;
   for (i = 3; i <= n; i++)
      if (prime[i])  {
         nprimes++;
      }
   printf("The number of primes found was %d\n",nprimes);
   printf("The total bases done was %d\n",totalwork);
}
