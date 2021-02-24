# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HEATED2D.

  Discussion:

    This code solves the steady state heat equation on a rectangular region.

    The sequential version of this program needs approximately
    18/epsilon iterations to complete. 


    The physical region, and the boundary conditions, are suggested
    by this diagram;

                   W = 0
             +------------------+
             |                  |
    W = 100  |                  | W = 100
             |                  |
             +------------------+
                   W = 100

    The region is covered with a grid of M by N nodes, and an N by N
    array W is used to record the temperature.  The correspondence between
    array indices and locations in the region is suggested by giving the
    indices of the four corners:

                  I = 0
          [0][0]-------------[0][N-1]
             |                  |
      J = 0  |                  |  J = N-1
             |                  |
        [M-1][0]-----------[M-1][N-1]
                  I = M-1

    The steady state solution to the discrete heat equation satisfies the
    following condition at an interior grid point:

      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )

    where "Central" is the index of the grid point, "North" is the index
    of its immediate neighbor to the "north", and so on.
   
    Given an approximate solution of the steady state heat equation, a
    "better" solution is given by replacing each interior point by the
    average of its 4 neighbors - in other words, by using the condition
    as an ASSIGNMENT statement:

      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )

    If this process is repeated often enough, the difference between successive 
    estimates of the solution will go to zero.

    This program carries out such an iteration, using a tolerance specified by
    the user, and writes the final estimate of the solution to a file that can
    be used for graphic processing.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 February 2021

  Author:

    Original C version by Michael Quinn.
    Modified C version by John Burkardt and Shirley Moore.

  Reference:

    Michael Quinn,
    Parallel Programming in C with MPI and OpenMP,
    McGraw-Hill, 2004,
    ISBN13: 978-0071232654,
    LC: QA76.73.C15.Q55.

  Parameters:

    double DIFF, the norm of the change in the solution from one iteration
    to the next.

    double MEAN, the average of the boundary values, used to initialize
    the values of the solution in the interior.

    double U[M][N], the solution at the previous iteration.

    double W[M][N], the solution computed at the latest iteration.
*/
{
# define M 500
# define N 500

  double diff;
  double epsilon = 0.001;
  int i;
  int iterations;
  int iterations_print;
  int j;
  double mean;
  double u[M][N];
  double w[M][N];
  double wtime;

  printf ( "\n" );
  printf ( "HEAT2D\n" );
  printf ( "  Sequential C version\n" );
  printf ( "  A program to solve for the steady state temperature distribution\n" );
  printf ( "  over a rectangular plate.\n" );
  printf ( "\n" );
  printf ( "  Spatial grid of %d by %d points.\n", M, N );
  printf ( "  The iteration will be repeated until the change is <= %e\n", epsilon ); 
/*
  Set the boundary values, which don't change. 
*/
  mean = 0.0;

    for ( i = 1; i < M - 1; i++ )
    {
      w[i][0] = 100.0;
    }
    for ( i = 1; i < M - 1; i++ )
    {
      w[i][N-1] = 100.0;
    }
    for ( j = 0; j < N; j++ )
    {
      w[M-1][j] = 100.0;
    }
    for ( j = 0; j < N; j++ )
    {
      w[0][j] = 0.0;
    }
/*
  Average the boundary values, to come up with a reasonable
  initial value for the interior.
*/
    for ( i = 1; i < M - 1; i++ )
    {
      mean = mean + w[i][0] + w[i][N-1];
    }
    for ( j = 0; j < N; j++ )
    {
      mean = mean + w[M-1][j] + w[0][j];
    }
  mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
  printf ( "\n" );
  printf ( "  MEAN = %f\n", mean );
/* 
  Initialize the interior solution to the mean value.
*/
    for ( i = 1; i < M - 1; i++ )
    {
      for ( j = 1; j < N - 1; j++ )
      {
        w[i][j] = mean;
      }
    }
/*
  iterate until the  new solution W differs from the old solution U
  by no more than EPSILON.
*/
  iterations = 0;
  iterations_print = 1;
  printf ( "\n" );
  printf ( " Iteration  Change\n" );
  printf ( "\n" );
  wtime = omp_get_wtime ( );

  diff = epsilon;

  while ( epsilon <= diff )
  {
/*
  Save the old solution in U.
*/
      for ( i = 0; i < M; i++ ) 
      {
        for ( j = 0; j < N; j++ )
        {
          u[i][j] = w[i][j];
        }
      }
/*
  Determine the new estimate of the solution at the interior points.
  The new solution W is the average of north, south, east and west neighbors.
*/
      for ( i = 1; i < M - 1; i++ )
      {
        for ( j = 1; j < N - 1; j++ )
        {
          w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;
        }
      }
    diff = 0.0;
      for ( i = 1; i < M - 1; i++ )
      {
        for ( j = 1; j < N - 1; j++ )
        {
          if ( diff < fabs ( w[i][j] - u[i][j] ) )
          {
            diff = fabs ( w[i][j] - u[i][j] );
          }
        }
      }

    iterations++;
    if ( iterations == iterations_print )
    {
      printf ( "  %8d  %f\n", iterations, diff );
      iterations_print = 2 * iterations_print;
    }
  } 
  wtime = omp_get_wtime ( ) - wtime;

  printf ( "\n" );
  printf ( "  %8d  %f\n", iterations, diff );
  printf ( "\n" );
  printf ( "  Error tolerance achieved.\n" );
  printf ( "  Wallclock time = %f seconds\n", wtime );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HEAT2D:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;

# undef M
# undef N
}
