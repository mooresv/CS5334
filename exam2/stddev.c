/*
 * stddev.c
 * Computes the standard deviation of an array of elements in parallel using MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

// Creates an array of random numbers. Each number has a value between 0 and 1
float *create_rand_nums(int num_elements) {
  float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (float)RAND_MAX);
  }
  return rand_nums;
}

int main(int argc, char** argv) {
  int rank, nprocs;
  int i;
  float local_sum, global_sum;
  float lsum_sq_diff, sum_sq_diff;
  float *rand_nums = NULL;

  if (argc != 2) {
    fprintf(stderr, "Usage: stddev num_elements_per_proc\n");
    exit(1);
  }

  int num_elements_per_proc = atoi(argv[1]);

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Create a random array of elements on each process
  srand(time(NULL)*rank);
  rand_nums = create_rand_nums(num_elements_per_proc);

  // Sum the numbers locally
  local_sum = 0;
  for (i = 0; i < num_elements_per_proc; i++) {
    local_sum += rand_nums[i];
  }

  // Reduce local sums into the global sum, leaving the result on every process

  /* Insert MPI collective call here */
   
  // Every process computes the mean
  float mean = global_sum / (num_elements_per_proc * nprocs);

  // Compute the local sum of the squared differences from the mean
  lsum_sq_diff = 0;
  for (i = 0; i < num_elements_per_proc; i++) {
    lsum_sq_diff += (rand_nums[i] - mean) * (rand_nums[i] - mean);
  }

  // Reduce the global sum of the squared differences to the root process

  /* Insert MPI collective call here */

  // The standard deviation is the square root of the mean of the squared
  // differences.
  if (rank == 0) {
    float stddev = sqrt(sum_sq_diff /
                        (num_elements_per_proc * nprocs));
    printf("Mean = %f, Standard deviation = %f\n", mean, stddev);
  }

  free(rand_nums);

  MPI_Finalize();
}
