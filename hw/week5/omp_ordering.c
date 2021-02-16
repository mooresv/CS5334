/*
 * OpenMP memory reordering example
*/

#include <stdio.h>
#include <omp.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

int main() {
    int a, b, c, d, i;
    int iterations = 0;
    int reorderings = 0;

#pragma omp parallel
{
    unsigned short xi[3];   // PRNG state variable
    #pragma omp for 
    for (i = 0; i < 2; i++) {
	// Read from /dev/urandom
        int fd = open("/dev/urandom", O_RDONLY);
        if (fd < 0) {
          perror("open /dev/urandom");
        exit(EXIT_FAILURE);
        }
        if (read(fd, xi, sizeof xi) != sizeof xi) {
            perror("read");
            exit(EXIT_FAILURE);
        }
        close(fd);
    }

    while (1) {
	c = d = 0;
	#pragma omp sections
	{
	#pragma omp section
	    { sleep(erand48(xi)); a = 1; c = b; }
	#pragma omp section
	    { sleep(erand48(xi)); b = 1; d = a; }
	}
	#pragma omp single
	{
	    iterations++;
            if (c == 0 && d == 0) {
	        reorderings++;
	        printf("%d reorderings in %d interations\n", reorderings, iterations);
	    }
	}
    }
}
}



