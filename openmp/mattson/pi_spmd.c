#include <omp.h>
#include <stdio.h>
static long num_steps = 1000000000; 
double step;
#define MAX_THREADS 8

void main () {
    int i, nthreads; double pi, sum[MAX_THREADS];
    step = 1.0/(double) num_steps;
    double start, end;

    start = omp_get_wtime();
    #pragma omp parallel
    {
        int i, id,nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if (id == 0) nthreads = nthrds;
        for (i=id, sum[id]=0.0;i< num_steps; i=i+nthrds) {
            x = (i+0.5)*step;
            sum[id] += 4.0/(1.0+x*x);
        }
     }
     end = omp_get_wtime();
     printf("Time for parallel region: %f\n", end-start);
     for(i=0, pi=0.0;i<nthreads;i++)pi += sum[i] * step;
}
