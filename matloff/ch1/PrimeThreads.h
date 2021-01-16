#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>  // required for threads usage

#define MAX_N 1000000000
#define MAX_THREADS 28

// shared variables
int nthreads,  // number of threads (not counting main())
    n,  // range to check for primeness
    prime[MAX_N+1],  // in the end, prime[i] = 1 if i prime, else 0
    nextbase;  // next sieve multiplier to be used
    // lock for the shared variable nextbase
    pthread_mutex_t nextbaselock = PTHREAD_MUTEX_INITIALIZER;
    // ID structs for the threads
//    pthread_t id[MAX_THREADS];
