/* File:    trap.c
 * Purpose: Calculate definite integral using trapezoidal 
 *          rule.
 *
 * Input:   a, b, n
 * Output:  Estimate of integral from a to b of f(x)
 *          using n trapezoids.
 *
 * Compile: gcc -g -Wall -o trap trap.c
 * Usage:   ./trap
 *
 * Note:    The function f(x) is hardwired.
 *
 * IPP:     Section 3.2.1 (pp. 94 and ff.) and 5.2 (p. 216)
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>

double f(double x);    /* Function we're integrating */
double Trap(double a, double b, int n, double h);

int main(void) {
   double  integral;   /* Store result in integral   */
   double  a, b;       /* Left and right endpoints   */
   int     n;          /* Number of trapezoids       */
   double  h;          /* Height of trapezoids       */
   clock_t start, end;

   printf("Enter a, b, and n, one per line\n");
   scanf("%lf", &a);
   scanf("%lf", &b);
   scanf("%d", &n);

   h = (b-a)/n;
   start = clock();
   integral = Trap(a, b, n, h);
   end = clock();
   printf("Elapsed time for Trap: %f seconds\n", ((double) (end-start))/CLOCKS_PER_SEC);
   
   printf("With n = %d trapezoids, our estimate\n", n);
   printf("of the integral from %f to %f = %.15f\n",
      a, b, integral);

   return 0;
}  /* main */

/*------------------------------------------------------------------
 * Function:    Trap
 * Purpose:     Estimate integral from a to b of f using trap rule and
 *              n trapezoids
 * Input args:  a, b, n, h
 * Return val:  Estimate of the integral 
 */
double Trap(double a, double b, int n, double h) {
   double integral;
   int k;

   integral = (f(a) + f(b))/2.0;
   for (k = 1; k <= n-1; k++) {
     integral += f(a+k*h);
   }
   integral = integral*h;

   return integral;
}  /* Trap */

/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) {
   double return_val;

   return_val = sqrt(fabs(sin(x)*sin(x)));
   return return_val;
}  /* f */
