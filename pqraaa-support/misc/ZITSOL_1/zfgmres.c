#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "./LIB/zheads.h"
#include "./LIB/zprotos.h"
#include "./LIB/zdefs.h"

#define  epsmac  1.0e-16


int zfgmres(SMatptr Amat, SPreptr PreMat, complex double *rhs, 
	    complex double *sol, double tol, int im, int *itmax, FILE *fp){ 
/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***                  
+-----------------------------------------------------------------------
| This is a simple Complex version of the ARMS preconditioned FGMRES algorithm. 
+-----------------------------------------------------------------------
| on entry:
|========== 
|
|(amat)   = matrix A in SpaFmt format. 
|
|(PreMat) = arms preconditioner structure consisting of:
|      (levmat) = precon lev. struture  containing linked list of 
|       LU(B), E, F matrices -- see LIB/heads.h for details/ 
|      (ilusch) = LU factorization of Schur complement matrix
|                                                                      
| rhs     = Complex vector of length n containing the right hand side.
|           
| sol     = Complex vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| im      = Krylov subspace dimension 
| itmax   = max number of iterations allowed. 
|
| on return:
|==========
| fgmresC  int =  0 --> successful return.
|          int =  1 --> convergence not achieved in itmax iterations.
|
| sol   == contains an approximate solution (upon successful return).
| itmax == has changed. It now contains the number of steps required
|          to converge -- 
+-----------------------------------------------------------------------
| work arrays: allocated internally and dynamically 
|=============
| vv    == work array of length [im+1][n] (used to store the Arnoldi
|          basis)
| hh    == work array of length [im][im+1] (Householder matrix)
| z     == work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
| armsol2, lusolD - preconditionning operation 
| BLAS1  routines.
|
+---------------------------------------------------------------------*/
  int n=Amat->n, maxits = *itmax; 
  int i, i1, ii, j, k, k1, its, retval, tmp1 = 1;
  complex double **hh,  *s, *rs, t1, **vv, **z, negt, rot;
  double t, beta, eps1=0.0, *c; 

   its = 0;
   vv = (complex double **) Malloc((im+1)*sizeof(complex double *), "fgmresC:1" );
   z  = (complex double **) Malloc(im*sizeof(complex double *), "fgmresC:2" );
   hh = (complex double **) Malloc(im*sizeof(complex double *), "fgmresC:3" );
   for (i=0; i<im; i++) {
     vv[i] = NULL;
     z[i]  = NULL;
   }
   vv[im] = NULL;
   for (i=0; i<im; i++) {
     hh[i] = (complex double *) Malloc((i+2)*sizeof(complex double), "fgmresC:4" );
   }
   c = (double *) Malloc(im*sizeof(double ), "fgmresC:5" );
   s = (complex double *) Malloc(im*sizeof(complex double), "fgmresC:6" );
   rs = (complex double *) Malloc((im+1)*sizeof(complex double), "fgmresC:7" );
   
   
/*-------------------------------------------------------------
|   outer loop starts here
+------------------------------------------------------------*/
label20:
/*-------------------------------------------------------------
|   compute initial residual vector
+------------------------------------------------------------*/
   if( vv[0] == NULL )
       vv[0] = (complex double *)Malloc( n*sizeof(complex double), "fgmresC:8" );
   Amat->zmatvec(Amat, sol, vv[0]);
   for (j=0; j<n; j++)
     vv[0][j] = rhs[j] - vv[0][j];
/*-------------------------------------------------------------
|  now vv[0] is the initial residual vector
+------------------------------------------------------------*/
   beta = dznrm2(&n, vv[0], &tmp1);
   /*   print resid. info if fp not null */
   if (fp != NULL) 
     fprintf(fp,"      %d       %e\n",its,beta) ;     
   if (beta == 0.0) goto label990;
   t = 1.0 / beta;
/*----------------------------------
|     normalize:  vv[0] = vv[0] / beta
+---------------------------------*/
   for (j=0; j<n; j++)
     vv[0][j] = vv[0][j]*t;
   if (its == 0) eps1 = tol*beta;
/*    ** initialize 1-st term  of rhs of hessenberg system   */
   rs[0] = beta + 0.0*I;
   i = -1;
 label4:
   i++;
   its++;
   i1 = i + 1;
/*------------------------------------------------------------
|    PRECONDITIONING    z_{j} = M^{-1} v_{j}
|                           w = A z_{j} = A M^{-1} v_{j}
+-----------------------------------------------------------*/
   if( z[i] == NULL )
     z[i] = (complex double *)Malloc( n*sizeof(complex double), "fgmresC:9" );
   
   if(PreMat == NULL)
   	memcpy(z[i],vv[i],n*sizeof(complex double));
   else{
/*-------------------- ## debug for nopre just uncomment next line! */
   		PreMat->zprecon(vv[i], z[i], PreMat) ;
	}

   if( vv[i1] == NULL) 
     vv[i1] = (complex double *)Malloc( n*sizeof(complex double), "fgmresC:10" );
   Amat->zmatvec(Amat, z[i], vv[i1]); 
/*------------------------------------------------------------
|     modified gram - schmidt...
|     h_{i,j} = (w,v_{i})
|     w  = w - h_{i,j} v_{i}
+------------------------------------------------------------*/
   for (j=0; j<=i; j++) {
      t1 = zdotc(&n, vv[j], &tmp1, vv[i1], &tmp1);
      hh[i][j] = t1;
      negt = -t1;
      zaxpy(&n, &negt, vv[j], &tmp1, vv[i1], &tmp1);
   }
/*----------------------------------
|     h_{j+1,j} = ||w||_{2}
+---------------------------------*/
   t = dznrm2(&n, vv[i1], &tmp1);
   hh[i][i1] = t + 0.0*I;
   if (t == 0.0) goto label58;
   t = 1.0/t;
/*----------------------------------
|     v_{j+1} = w / h_{j+1,j}
+---------------------------------*/
   for (k=0; k<n; k++)
     vv[i1][k] = vv[i1][k]*t;
/*---------------------------------------------------
|     done with modified gram schimdt and arnoldi step
|     now  update factorization of hh
+--------------------------------------------------*/
 label58:
/*--------------------------------------------------------
|   perform previous transformations  on i-th column of h
+-------------------------------------------------------*/
   for (k=1; k<=i; k++) {
     k1 = k-1;
     t1 = hh[i][k1];
     hh[i][k1] = (c[k1])*t1 + s[k1]*hh[i][k];
     hh[i][k] = -conj(s[k1])*t1 + (c[k1])*hh[i][k];
    }
      
/*---------------------------------------------------
|     if gamma is zero then any small value will do...
|     will affect only residual estimate
+--------------------------------------------------*/
   //  if (gam == 0.0) gam = epsmac;
/*---------------------------------------------------
|     get  next plane rotation
+--------------------------------------------------*/

  zclartg(hh[i][i], hh[i][i1], &c[i], &s[i], &rot);
  // printf("c1 = %f \n", c[i]);
 
     rs[i1] = -conj(s[i])*rs[i];
     rs[i] =  (c[i])*rs[i];
/*----------------------------------------------------
|   determine residual norm and test for convergence
+---------------------------------------------------*/
//   hh[i][i] = conj(c[i])*hh[i][i] + s[i]*hh[i][i1];
	hh[i][i] = rot;
   beta = cabs(rs[i1]);
   if (fp != NULL)
   fprintf(fp,"      %d       %e\n",its,beta) ;     
   if ( (i < im-1) && (beta > eps1) && (its < maxits) )  goto label4;
/*---------------------------------------------------
|     now compute solution
|     first solve upper triangular system
+--------------------------------------------------*/
   rs[i] = rs[i]/hh[i][i];
   for (ii=1; ii<=i; ii++) {
     k=i-ii;
     k1 = k+1;
     t1=rs[k];
     for (j=k1; j<=i; j++){
       t1 = t1 - hh[j][k]*rs[j];
     	 rs[k] = t1/hh[k][k];
     	 }
   }
/*-----------------------------------------------------
|   form linear combination of v[i]'s to get solution   
+----------------------------------------------------*/
   for (j=0; j<=i; j++) {
     t1 = rs[j];
     for (k=0; k<n; k++)
       sol[k] += t1*z[j][k];
   }
/*----------------------------------------------------
|     restart outer loop  when enecessary  
+---------------------------------------------------*/
   if (beta <= eps1) goto label990;
   if (its >= maxits) goto label991;
   goto label20;
label990:
   retval = 0;
   goto label888;
label991:
   retval = 1;
label888:
   for (i=0; i<=im; i++)
     if( vv[i] ) free(vv[i]);
   free(vv);
   for (i=0; i<im; i++){
     free(hh[i]);
     if( z[i] ) free(z[i]);
   }
   free(hh);
   free(z);
   free(c);
   free(s);
   free(rs);
   *itmax = its; 
   return retval;
}
/*-----------------end of fgmresC ------------------------------------
----------------------------------------------------------------------*/
