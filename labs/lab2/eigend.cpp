/* Copyright (c) 2013, Intel Corporation
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, 
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, 
 *   this list of conditions and the following disclaimer in the documentation 
 *   and/or other materials provided with the distribution.
 * - Neither the name of Intel Corporation nor the names of its contributors 
 *   may be used to endorse or promote products derived from this software 
 *   without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <math.h>

int eigend(double aptr[],double rptr[],int nsize)
/*
 *  EIGENSOLVER FOR SMALL MATRICES
 *  Solve for eigenvalues and eigenvectors of a[][]
 *  using threshold Jacobi technique.
 *  a[][] must be a real floating point symmetric matrix.
 *  Input matrix a[][] is destroyed in the process
 *  and the diagonal elements of a[][] are loaded
 *  with its eigenvalues.  The matrix r[][] is 
 *  the output matrix of eigenvectors.
 *  The input integer nsize determines size of matrices a & r.
 *  The number of annihilations performed is returned.
 *  Paul J. Besl   5-30-84
 */
{
double anorm,wnorm,vnorm,epsilon,final,threshold;
double sin0,cos0,tan0,cot20,h,tau;
double aii,ajj,aij,aki,akj,rki,rkj;
double *a,*r;
int nsqr,count;
int i,j,k,n;
/*--------------------------------------------------------*/
   n = nsize; nsqr = n*n; a = aptr ; r = rptr;
   vnorm = 0.0; wnorm = 0.0;
   epsilon = 1.0e-24;  count=0;
/*
 * initialize r[][] as an nxn identity matrix
 * compute on- and off-diagonal norms of a[][]
 */
   for(k=0; k<nsqr; ++k)
   { /* row index=i and column index=j */
     i=k/n;  j=k-n*i;
     if(i == j)  /* on-diagonal */
     { *r++ = 1.0; aii = *a++; wnorm  += aii*aii; }
     else  /* off-diagonal */
     { *r++ = 0.0; aij = *a++; vnorm  += aij*aij; }
   }
/*
 * compute matrix norms
 */
   anorm = vnorm + wnorm;  anorm = sqrt(anorm);
   vnorm = sqrt(vnorm);    wnorm = sqrt(wnorm);
/*
 * compute initial threshold and final threshold
 */
   threshold = (vnorm) / ((float) n);
   final = epsilon * threshold;
   if( final < 1.0e-31 ) final = 1.0e-31;
   a = aptr ; r = rptr;

/* main computation loop:  decrease threshold */ 
/* ========================================== */
   while ( threshold > final )
   {  threshold = threshold / ((float) n);

   /* loop over upper off-diagonal elements */
   /* row index=i and column index=j */
      for (i=0;  i<n-1; ++i)
      { for(j=i+1; j<n; ++j)
        { aij = *(a + j + i*n);
          if(fabs(aij) > threshold)
          {
          /* compute sine, cosine, tangent, & cotangent */
            aii = *(a + i + i*n );
            ajj = *(a + j + j*n );
            cot20 = (ajj - aii)/(2.0*aij);
            if( fabs(cot20) < 1.0e+15 )
            { tan0  = fabs( cot20 ) + sqrt( 1. + cot20*cot20 ); }
            else { tan0 = 2.0*fabs( cot20 ); }
            tan0  = 1.0/tan0;
            if(cot20 < 0) tan0 = -tan0;
            cos0  = 1.0/sqrt(1.0 + tan0*tan0);
            sin0  = cos0*tan0;
            tau   = sin0/(1.0 + cos0);
            h     = tan0 * aij;
          /*
           * update two diagonal and two zeroed-out terms
           */
            *(a + i + i*n) = aii - h;
            *(a + j + j*n) = ajj + h;
            *(a + i + j*n) = 0.0;
            *(a + j + i*n) = 0.0;
            ++count;
          /*
           * update "a" in two columns and two rows
           * update "r" in two columns
           * k = row index of updated columns
           */
            for(k=0; k<n; ++k)
            { /* always rotate r matrix */
              rki= *(r + i + k*n);
              rkj= *(r + j + k*n);
 	          *(r + i + k*n) += (-sin0)*(rkj+(tau*rki));
              *(r + j + k*n) += ( sin0)*(rki-(tau*rkj));
              if(k != i && k != j)
              { /* update a matrix unless k=i or k=j */
                akj= *(a + j + k*n);
                aki= *(a + i + k*n);
 		        *(a + i + k*n) += (-sin0)*(akj+(tau*aki));
                *(a + j + k*n) += ( sin0)*(aki-(tau*akj));
 		        *(a + k + i*n) = *(a + i + k*n); 
                *(a + k + j*n) = *(a + j + k*n);
              } /* end of if k=i or j */
            }  /* end of k-loop update */
          }   /* end of if aij > threshold */
        }    /* end of j-loop */
      }     /* end of i-loop */
   }       /* end of threshold loop */
/*
 * round-off small numbers below final to zero
 */
#if 1
   a = aptr ; r = rptr;
   for(k=0;k<nsqr;++k)
   { if( fabs( *a ) < final )  *a=0.0; 
     if( fabs( *r ) < final )  *r=0.0;
     ++a; ++r;
   }
#endif
   return(count);
/*----------------------------------------------------*/
}

