#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "./zheads.h"
#include "./zprotos.h"
#include "./zdefs.h"


/*---------------------------------------------------------------------------*
 * Complex version of routines for performing matrix operations - MatOps
 *---------------------------------------------------------------------------*/

void *Malloc( int, char * );

void zmatvec(csptr mata, complex double *x, complex double *y)
{
/*---------------------------------------------------------------------
| This function does the matrix vector product y = A x.
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| x     = a vector
|
| on return
| y     = the product A * x
|--------------------------------------------------------------------*/
/*-------------------- local variables                               */
  int i, k, *ki;
  complex double *kr;
  for (i=0; i<mata->n; i++) {
    y[i] = 0.0 + 0.0*I;
    kr = mata->ma[i];
    ki = mata->ja[i];
    for (k=0; k<mata->nzcount[i]; k++) 
      y[i] += kr[k] * x[ki[k]];
	}
  return;
}
/*---------------end of matvec------------------------------------------
----------------------------------------------------------------------*/

void zLsol(csptr mata, complex double *b, complex double *x)
{
/*---------------------------------------------------------------------
| This function does the forward solve L x = b.
| Can be done in place.
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| b     = a vector
|
| on return
| x     = the solution of L x = b 
|--------------------------------------------------------------------*/
/*   local variables    */
  int i, k;
  complex double *kr;
  int *ki;
  for (i=0; i<mata->n; i++) {
    x[i] = b[i];
    if ( mata->nzcount[i] > 0 ) {
      kr = mata->ma[i];
      ki = mata->ja[i];
      for (k=0; k<mata->nzcount[i]; k++)
	x[i] -= kr[k]*x[ki[k]];
    }
  }
  return;
}
/*---------------end of Lsol-----------------------------------------
----------------------------------------------------------------------*/
void zUsol(csptr mata, complex double *b, complex double *x)
{
/*---------------------------------------------------------------------
| This function does the backward solve U x = b.
| Can be done in place.
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| b    = a vector
|
| on return
| x     = the solution of U * x = b 
|
|---------------------------------------------------------------------*/
/*   local variables    */
  int i, k, *ki;
  complex double *kr;
  for (i=mata->n-1; i>=0; i--) {
    kr = mata->ma[i];
    ki = mata->ja[i];
    x[i] = b[i] ;
    for (k=1; k<mata->nzcount[i]; k++)
      x[i] -= kr[k] * x[ki[k]];
    x[i] *= kr[0];
  }
  return;
}
/*----------------end of Usol----------------------------------------
----------------------------------------------------------------------*/
int zdescend(p4ptr levmat, complex double *x, complex double *wk)
{
/*---------------------------------------------------------------------
| This function does the (block) forward elimination in ARMS
|                       new       old
|     |            |  |     |    |    |
|     | L        0 |  | wx1  |    | x1 |
|     |            |  |     | =  |    | 
|     | EU^{-1}  I |  | wx2  |    | x2 |
|     |            |  |     |    |    |
| x used and not touched -- or can be the same as wk.
|--------------------------------------------------------------------*/
/*  local variables   */
  int j, len=levmat->n, lenB=levmat->nB, *iperm=levmat->rperm; 
  complex double *work = levmat->wk; 
/*------------------------------------------------------
|   apply permutation P to rhs 
|-----------------------------------------------------*/
  for (j=0; j<len; j++)
    work[iperm[j]] = x[j] ;
  zLsol(levmat->L, work, wk);       /* sol:   L x = x                 */
  zUsol(levmat->U, wk, work);        /* sol:   U work(2) = work         */
/*-------------------- compute x[lenb:.] = x [lenb:.] - E * work(1) */
  zmatvecz (levmat->E, work, &work[lenB], &wk[lenB]) ; 
  return 0;
}
/*----end-of-descend---------------------------------------------------
|----------------------------------------------------------------------
|--------------------------------------------------------------------*/
int zascend (p4ptr levmat, complex double *x, complex double *wk) 
{
/*---------------------------------------------------------------------
| This function does the (block) backward substitution: 
|
|     |            |  |     |    |    |
|     | U  L^{-1}F |  | wk1 |    | x1 |
|     |            |  |     | =  |    |
|     | 0       S  |  | wk2 |    | x2 |  <<-- x2 already computed.
|     |            |  |     |    |    |       and we need x1
|
|    with x2 = S^{-1} wk2 [assumed to have been computed ] 
|--------------------------------------------------------------------*/
 /*--------------------  local variables  */
  int j, len=levmat->n, lenB=levmat->nB, *qperm=levmat->perm;
  complex double *work = levmat->wk; 
  /*-------------------- copy x onto wk */  
  zmatvec(levmat->F, &x[lenB], work);   /*  work = F * x_2   */
  zLsol(levmat->L, work, work);         /*  work = L \ work    */
  for (j=0; j<lenB; j++)               /*  wk1 = wk1 - work  */
    work[j] = x[j] - work[j];
  zUsol(levmat->U, work, work);         /*  wk1 = U \ wk1 */ 
  memcpy(&work[lenB],&x[lenB],(len-lenB)*sizeof(complex double));
/*---------------------------------------
|   apply reverse permutation
|--------------------------------------*/
  for (j=0; j<len; j++)
     wk[j] = work[qperm[j]];     
  return 0;
}
/*----end-of-ascend----------------------------------------------------
|----------------------------------------------------------------------
|--------------------------------------------------------------------*/

void zmatvecz(csptr mata, complex double *x, complex double *y, complex double *z) 
{
/*---------------------------------------------------------------------
| This function does the matrix vector  z = y - A x.
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| x, y   = two input vector
|
| on return
| z    = the result:  y - A * x
| z-location must be different from that of x 
| i.e., y and x are used but not modified.
|--------------------------------------------------------------------*/
/*   local variables    */
  int i, k, *ki;
  complex double *kr, t;
  for (i=0; i<mata->n; i++) {
    kr = mata->ma[i];
    ki = mata->ja[i];
    t = y[i] ;
    for (k=0; k<mata->nzcount[i]; k++)
      t -= kr[k] * x[ki[k]];
    z[i] = t; 
  }
  return;
}
/*---------------end of matvecz----------------------------------------
 *--------------------------------------------------------------------*/

p4ptr zLvsol2(complex double *x, int nlev, p4ptr levmat, ilutptr ilusch) 
{
  /* Macro L-solve -- corresponds to left (L) part of arms
  |  preconditioning operation -- 
  |  on entry : 
  |   x =  right- hand side to be operated on by the preconditioner
  |  on return : x is overwritten
  |   x =  output result of operation 
  |  
  |  Note : in-place operation -- b and x can occupy the same space..
  | --------------------------------------------------------------------*/ 
/*-------------------- local variables  */
  int nloc=levmat->n, first, lenB;
  p4ptr last=levmat; 
/*-------------------- take care of  special cases :  nlev==0 --> lusol  */
  if (nlev == 0) {
    zSchLsol(ilusch,x);
    return (last);
  }
  first = 0;
/*-------------------- descend                                      */
  while (levmat) { 
    nloc=levmat->n;
    lenB =levmat->nB;
/*-------------------- left scaling                                  */
    if (levmat->D1 !=  NULL) 
      zdscale(nloc,levmat->D1, &x[first],  &x[first]); 
/*--------------------  RESTRICTION/ DESCENT OPERATION  */
    if (lenB) 
      zdescend (levmat, &x[first], &x[first]);
    first += lenB; 
    last = levmat;
    levmat = levmat->next;
/*---------------------------------------------------------------------
| next level 
+--------------------------------------------------------------------*/
   } 
  zSchLsol(ilusch,&x[first]);
  return last; 
}

int zUvsol2(complex double *x, int nlev, int n, p4ptr levmat,
	   ilutptr ilusch) 
{
  /* Macro U-solve -- corresponds to left (L) part of arms
  |  preconditioning operation -- 
  |  on entry : 
  |  b  =  right- hand side to be operated on by the preconditioner
  |  on return  = x has been overwritten =
  |  x  =  output result of operation 
  |  
  |  Note : in-place operation -- b and x  can occupy the same space..
  | --------------------------------------------------------------------*/ 

/*-------------------- local variables  */
  int nloc, lenB, first; 
  /*-------------------- work array                                        */
  /*-------------------- take care of  special cases :  nlev==0 --> lusol  */
/*-------------------- case of zero levels                             */
  if (nlev == 0) { 
    zSchUsol(ilusch, x);  
    return(0);
  }
/*-------------------- general case                               */
  nloc=levmat->n; 
  lenB=levmat->nB; 
  first = n - nloc; 
/*-------------------- last level                                 */
  first += lenB; 
  zSchUsol(ilusch, &x[first]); 
/*-------------------- other levels                               */
  while (levmat) {
    nloc = levmat->n; 
    first -= levmat->nB;
    if (levmat->n) 
      zascend(levmat, &x[first],&x[first]);
/*-------------------- right scaling */
    if (levmat->D2 !=  NULL) 
      zdscale(nloc, levmat->D2, &x[first], &x[first]) ;
    levmat = levmat->prev; 
  }
  return 0;
/*--------------------  PROLONGATION/ ASCENT OPERATION */
}

int zarmsol2(complex double *x,  arms Prec) 
{ /* combined preconditioning operation -- combines the
  |  left and right actions. 
  | 
  |  on entry : 
  |   x =  right- hand side to be operated on by the preconditioner
  |  on return : x is overwritten - 
  |   x =  output result of operation 
  |  
  |  Note : in-place operation -- b and x can occupy the same space..
  | --------------------------------------------------------------------*/ 
/*-------------------- local variables  */
  p4ptr levmat = Prec->levmat;
  ilutptr ilusch = Prec->ilus;
  int nlev = Prec->nlev;
  int n = levmat->n; 
  p4ptr last;

  if (nlev == 0) {
    n = ilusch->n;
    zSchLsol(ilusch, x); 
    zSchUsol(ilusch, x); 
    return 0;
  }
  last = zLvsol2(x, nlev, levmat, ilusch) ;
  zUvsol2(x, nlev, n, last, ilusch) ; 
  return 0; 
}

void zSchLsol(ilutptr ilusch, complex double *y) 
{
/*---------------------------------------------------------------------
|  Forward solve for Schur complement part = 
|----------------------------------------------------------------------
| on entry:
| ilusch  = the LU matrix as provided from the ILU functions.
| y       = the right-hand-side vector
|
| on return
| y       = solution of LU x = y. [overwritten] 
|---------------------------------------------------------------------*/
/*-------------------- local variables                        */
  int n = ilusch->n, j, *perm = ilusch->rperm;
  complex double *work = ilusch->wk; 
/*-------------------- begin: right scaling                          */
   if (ilusch->D1 != NULL) 
     zdscale(n, ilusch->D1, y, y); 
/*-------------------- ONE SIDED ROW PERMS */
   if (perm != NULL) { 
     for (j=0; j<n; j++)
       work[perm[j]] = y[j]; 
/*--------------------  L solve proper */
     zLsol(ilusch->L, work, y); 
   } else 
     zLsol(ilusch->L, y, y); 
/*---------------end of SchLsol---------------------------------------
----------------------------------------------------------------------*/
}

void zSchUsol(ilutptr ilusch, complex double *y) 
{
/*---------------------------------------------------------------------
| U-solve for Schur complement  - 
|----------------------------------------------------------------------
| on entry:
| ilusch  = the LU matrix as provided from the ILU functions.
| y       = the right-hand-side vector
|
| on return 
| y       = solution of U x = y. [overwritten on y] 
|----------------------------------------------------------------------*/
  int n = ilusch->n, j,  *perm = ilusch->perm, *cperm;
  complex double *work = ilusch->wk; 
/* -------------------- begin by U-solving */
/*-------------------- CASE: column pivoting  used (as in ILUTP) */
  if (ilusch->perm2 != NULL) {
    zUsol(ilusch->U, y, y);
    cperm = ilusch->perm2; 
    for (j=0; j<n; j++)
      work[cperm[j]] = y[j];
  }
  else
/*-------------------- CASE: no column pivoting  used                   */
    zUsol(ilusch->U, y, work);
/*-------------------- generic permutation                              */
  if (perm != NULL) {
    for (j=0; j<n; j++)
      y[j] = work[perm[j]];
  } else
    memcpy(y, work,n*sizeof(complex double));
    
/*-------------------- case when diagonal scaling is done on columns    */ 
    if (ilusch->D2 !=NULL) 
      zdscale(n, ilusch->D2, y, y);
}
/*---------------end of SchUsol---------------------------------------
----------------------------------------------------------------------*/

int zlusolC( complex double *y, complex double *x, iluptr lu ) 
{
/*----------------------------------------------------------------------
 *    performs a forward followed by a backward solve
 *    for LU matrix as produced by iluk
 *    y  = right-hand-side 
 *    x  = solution on return 
 *    lu = LU matrix as produced by iluk. 
 *--------------------------------------------------------------------*/
    int n = lu->n, i, j, nzcount, *ja;
    complex double *D;
    csptr L, U;

    L = lu->L; 
    U = lu->U;
    D = lu->D;

    /* Block L solve */
    for( i = 0; i < n; i++ ) {
        x[i] = y[i];
        nzcount = L->nzcount[i];
        ja = L->ja[i];
        for( j = 0; j < nzcount; j++ ) {
            x[i] -= x[ja[j]] * L->ma[i][j];
        }
    }
    /* Block -- U solve */
    for( i = n-1; i >= 0; i-- ) {
        nzcount = U->nzcount[i];
        ja = U->ja[i];
        for( j = 0; j < nzcount; j++ ) {
            x[i] -= x[ja[j]] * U->ma[i][j];
        }
        x[i] *= D[i];
    }
    return (0); 
}

int zrpermC(csptr mat, int *perm)
{
/*----------------------------------------------------------------------
| This function permutes the rows of a matrix in SpaFmt format. 
| rperm  computes B = P A  where P is a permutation matrix.  
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination row number of row number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SpaFmt format.
|
|
| on return:
| ----------
| (amat) = P A stored in SpaFmt format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
  int **addj, *nnz, i, size=mat->n;
  complex double **addm;
  addj = (int **) Malloc(size*sizeof(int *), "rpermC:1" );
  addm = (complex double **) Malloc(size*sizeof(complex double *), "rpermC:2" );
  nnz = (int *) Malloc(size*sizeof(int), "rpermC:3" );
  for (i=0; i<size; i++) {
    addj[perm[i]] = mat->ja[i];
    addm[perm[i]] = mat->ma[i];
    nnz[perm[i]] = mat->nzcount[i];
  }
  for (i=0; i<size; i++) {
    mat->ja[i] = addj[i];
    mat->ma[i] = addm[i];
    mat->nzcount[i] = nnz[i];
  }
  free(addj);
  free(addm);
  free(nnz);
  return 0;
}
/*------- end of rperm ------------------------------------------------- 
|---------------------------------------------------------------------*/
int zcpermC(csptr mat, int *perm) 
{
/*----------------------------------------------------------------------
| This function permutes the columns of a matrix in SpaFmt format.
| cperm computes B = A P, where P is a permutation matrix.
| that maps column j into column perm(j), i.e., on return 
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination column number of column number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (mat) = a matrix stored in SpaFmt format.
|
|
| on return:
| ----------
| (mat) = A P stored in SpaFmt format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
   int i, j, *newj, size=mat->n, *aja;
   newj = (int *) Malloc(size*sizeof(int), "cpermC" );
   for (i=0; i<size; i++) {
      aja = mat->ja[i];
      for (j=0; j<mat->nzcount[i]; j++)
	 newj[j] = perm[aja[j]];
  
      for (j=0; j<mat->nzcount[i]; j++)
	 aja[j] = newj[j];
      mat->ja[i] = aja;
   }
   free(newj);
   return 0;
}
/*------- end of cperm ------------------------------------------------- 
|---------------------------------------------------------------------*/
int zdpermC(csptr mat, int *perm) 
{
/*----------------------------------------------------------------------
| This function permutes the rows and columns of a matrix in 
| SpaFmt format.  dperm computes B = P^T A P, where P is a permutation 
| matrix.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SpaFmt format.
|
|
| on return:
| ----------
| (amat) = P^T A P stored in SpaFmt format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
   if (zrpermC(mat, perm)) return 1;
   if (zcpermC(mat, perm)) return 1;
   return 0;
}
/*------- end of dperm ------------------------------------------------- 
|---------------------------------------------------------------------*/

int zcondestLU( iluptr lu, complex double *y, complex double *x, FILE *fp )
{
    int n = lu->n, i;
    double norm = 0.0;

    for( i = 0; i < n; i++ ) {
        y[i] = 1.0 + 0.0*I;  
    }
    zlusolC( y, x, lu );
    for( i = 0; i < n; i++ ) {
        norm = max( norm, cabs(x[i]) );
    }
    fprintf( fp, "ILU inf-norm lower bound : %16.2f\n", norm );
    if( norm > 1e30 ) {
        return -1;
    }
    return 0;
}

int zcondestArms(arms armspre, complex double *y, FILE *fp ){
/*-------------------- simple estimate of cond. number of precon */
  int n = armspre->n, i;
  double norm = 0.0;
  
  for( i = 0; i < n; i++ ) 
    y[i] = 1.0 + 0.0*I; 
  zarmsol2(y, armspre)  ;
  for( i = 0; i < n; i++ ) {
    norm = max( norm, cabs(y[i]) );
  }
  fprintf( fp, "ARMS inf-norm lower bound = : %16.2f\n", norm );
  if( norm > 1e30 ) {
    return -1;
  }
  return 0;
}

void zmatvecCSR(SMatptr mat, complex double *x, complex double *y){
  /*-------------------- matvec for csr format using the SMatptr struct*/
  zmatvec(mat->CSR, x, y )  ;
}

/*
iluc related 
void zmatvecLDU(SMatptr mat, complex double *x, complex double *y){
-------------------- matvec for arms format using the Smatrix struct
zlumatvec(mat->LDU, x, y )  ;
}
*/

/*-------------------- preconditioning operations */


int  zpreconILU(complex double *x, complex double *y, SPreptr mat){ 
  /*-------------------- precon for csr format using the SPre struct*/
  return zlusolC(x, y, mat->ILU)  ;
}


/*   
iluc related -- removed for now
int  zpreconLDU(complex double *x, complex double *y, SPreptr mat){
 -------------------- precon for LDU format using the SPre struct
 return zlumsolC(x, y, mat->ILU)  ;
}
*/

int  zpreconARMS(complex double *x, complex double *y, SPreptr mat){
  /*-------------------- precon for arms format using the SPre struct*/
  int n = (mat->ARMS)->n , i; 
  for(i = 0; i < n; i++ )
        y[i] = x[i];
 
  return zarmsol2(y, mat->ARMS)  ;
}

