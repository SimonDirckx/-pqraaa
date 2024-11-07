#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "./zheads.h"
#include "./zprotos.h"
#define  ALPHA  0.00001

/*---------------------------------------------------------------------------*
 * Complex version of PQ.c
 *---------------------------------------------------------------------------*/

int zPQperm(csptr mat, int bsize, int *Pord, int *Qord, int *nnod, 
	     double tol) {
/*--------------------------------------------------------------------- 
| Algorithm for nonsymmetric  block selection - 
|----------------------------------------------------------------------
|     Input parameters:
|     -----------------
|     (mat)  =  matrix in SpaFmt format
|     
|     tol    =  a tolerance for excluding a row from B block 
|
|     bsize not used here - it is used in arms2.. 
|
|     Output parameters:
|     ------------------ 
|     Pord   = row permutation array.  Row number i will become row number 
|		Pord[i] in permuted matrix. (Old to new labels) 
|     Qord   = column permutation array.  Column number j will become 
|		row number Qord[j] in permuted matrix.
|             [destination lists] - i.e., old to new arrays= 
|     
|     nnod   = number of elements in the B-block 
|     
|---------------------------------------------------------------------*/ 
/*--------------------   local variables   */
  int *icor, *jcor, *row; 
  int i, j, ii, k, col, jj, rnz, nzi, n=mat->n, count, numnode;
  double  aij, rn;
  complex double *mrow;  
/*-------------------- prototypes */
  int zpreSel(csptr mat, int *icor, int *jcor, int job, 
	     double tol, int *count); 
  int check_perm(int, int *) ;   /* debug only */
/*-----------------------------------------------------------------------*/  
  for (j=0; j<n; j++) {
    Pord[j] = -1;
    Qord[j] = -1; 
  }
  icor = (int *) malloc(n*sizeof(int));
  jcor = (int *) malloc(n*sizeof(int));
  if ( (icor==NULL) || (jcor==NULL) ) return 1;
  numnode = 0;
  count = 0;
/*-------------------- wDiag selects candidate entries in a sorted oder */
  i = 1; 
  zpreSel(mat, icor, jcor, i, tol, &count) ;

/* ---- DEBUG: Check selection --- */  
/*  int tt;
        printf("\n  icor	jcor\n"); 
        for(tt=0;tt<count; tt++)
        {
        		printf("	%d	%d\n", icor[tt], jcor[tt]);
        	}
        	printf("\n");
 */ 
        
/*-------------------- add entries one by one to diagonal */
  /* needs recoding so as to scan rows only once instead of 2 */
  for (i = 0; i<count; i++){ 
    ii = icor[i];
    jj = jcor[i]; 
    if (Qord[jj] != -1) continue;
/*-------------------- */
    row = mat->ja[ii]; 
    mrow = mat->ma[ii];
    nzi = mat->nzcount[ii] ;
/*-------------------- rnz = already assigned cols (either in L or F) */
    rn = cabs(mrow[0]); 
    rnz = (nzi-1) ; 
    for (k=0; k < nzi; k++) {
      aij = cabs(mrow[k]);
      col = row[k];
      if (Qord[col] >=0 ) {
	rn -= aij; 
	rnz-- ; 
      }
      else if (Qord[col] == -2) {
	rnz--;
      }
    } 
    if (rn < 0.0) continue;   
    Pord[ii] = numnode;
    Qord[jj] = numnode;
    numnode++; 
/*-------------------- acceptance test among others */    
    for (k=0; k < nzi; k++) {
      col = row[k];
      if (Qord[col] != -1) continue;
      aij = cabs(mrow[k]);
      if (rnz*aij > rn) 
	Qord[col] = -2;
      else 
	rn -= aij;
      rnz--;
    }
  }
  /*-------------------- number of B nodes */
  *nnod = numnode; 
  /* printf(" nnod found = %d \n",*nnod);  */ 
/*--------------------------------------------------
|    end-main-loop - complete permutation arrays
|-------------------------------------------------*/
  for (i=0; i<n; i++)
    if (Pord[i] < 0) 
      Pord[i] = numnode++;
  
  if (numnode != n) {
    printf("  ** counting error - type 1 \n"); return 1; }   
  numnode = *nnod;
  for (j=0; j<n; j++)
    if (Qord[j] < 0)
      Qord[j] = numnode++;
/*--------------------              */ 
  if (numnode != n) {
    printf(" **  counting error type 2 \n"); return 1; }

/*-------------------- debugging - check permutations */
  /* 
     printf(" checking P  and Q  :    -->  \n") ;
     check_perm(n, Pord) ;
     check_perm(n, Qord) ; 
  */
/*--------------------------------------------------
|  clean up before returning
|-------------------------------------------------*/
  free(icor);
  free(jcor);
  return 0;
}
/*---------------------------------------------------------------------
|-----end-of-zPQperm--------------------------------------------------
|--------------------------------------------------------------------*/

int zpreSel(csptr mat, int *icor, int *jcor, int job, double tol, int *count) 
{
/*---------------------------------------------------------------------
| does a preselection of possible diagonal entries. will return a list
| of "count" bi-indices representing "good" entries to be selected as 
| diagonal elements -- the order is important (from best to
| to worst). The list is in the form (icor(ii), jcor(ii)) 
|
|      ON ENTRY: 
|       mat   = matrix in csptr format 
|       tol   = tolerance used for selecting best rows -|
|       job   = indicates whether or not to permute the max entry in 
|               each row to first position 
|        NOTE: CAN RECODE BY HAVING JCOR CARRY THE INDEX IN ROW[I] WHERE
|              MAX IS LOCATED.. 
|       
|      ON RETURN: 
|       icor  = list of row indices of entries selected 
|       jcor  = list of column indices of entries selected 
|       count = number of entries selected (size of B block) 
|--------------------------------------------------------------------*/
  int i, k, kmax, n=mat->n, len, col, jmax, countL;
  int *nz, *jcol; 
  double rownorm,t, wmax,  tmax, *weight;
  complex double *mrow, t1;
  void qsortR2I(double *, int *i, int *, int, int);
/*--------------------begin */
/*-----------------------------------------------------------------------*/
  nz =mat->nzcount;
  len = 0; 
  weight = (double *) malloc(n*sizeof(double));
  if ( weight==NULL) return 1;  
  /*-------------------- compute max entry for each row */
  wmax = 0.0; 
  for (i=0; i<n; i++) {
    jcol = mat->ja[i];
    mrow = mat->ma[i];
    tmax = 0.0; 
    kmax = 0; 
    rownorm = 0.0; 
    for (k = 0; k<nz[i]; k++) 
    {
      col = jcol[k] ; 
      t = cabs(mrow[k]);  
      if (t != 0.0) 
      {
			rownorm += t; 
			if (tmax < t) 
			{
	  			tmax = t;
	  			kmax = k;
			}
      }
    } 
    jmax = jcol[kmax];
    jcor[i] = jmax; 
    if (job && kmax != 0) {
      t1 = mrow[kmax];
      mrow[kmax] = mrow[0];
      mrow[0] = t1;
      jcol[kmax] = jcol[0];
      jcol[0] = jmax;
    }
/*-------------------- save max diag. dominance ratio  */
    t = tmax / rownorm;  
    if (wmax < t)  wmax = t; 
    weight[i] = t;
    /* remove!! ALREADY ASSIGNED  */
    jcor[i] = jmax;
  }
/*-------------------- now select according to tol */
  countL = 0; 
  for (i=0; i<n; i++) {
    t = weight [i] ;
    col = jcor[i];
    if (t < wmax*tol) continue ;
    weight[countL] =  t /((double) nz[i]) ; 
    icor[countL] = i; 
    jcor[countL] = col;
    countL++;
  }
/*-------------------- sort them  */
  qsortR2I(weight, icor, jcor, 0, countL-1);
  *count = countL;
  
free(weight);
  return 0;
}
/*---------------------------------------------------------------------
|---- end of zpreSel ---------------------------------------------------
|--------------------------------------------------------------------*/
