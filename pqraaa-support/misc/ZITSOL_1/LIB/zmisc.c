#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "zheads.h"

/*---------------------------------------------------------------------------*
 * Complex version of skit
 *---------------------------------------------------------------------------*/

void *Malloc( int, char * );
int zqsplitC(complex double *a, int *ind, int n, int ncut)
{
/*----------------------------------------------------------------------
|     does a quick-sort split of a complex real array.
|     on input a[0 : (n-1)] is a real array
|     on output is permuted such that its elements satisfy:
|
|     abs(a[i]) >= abs(a[ncut-1]) for i < ncut-1 and
|     abs(a[i]) <= abs(a[ncut-1]) for i > ncut-1
|
|     ind[0 : (n-1)] is an integer array permuted in the same way as a.
|---------------------------------------------------------------------*/
   complex double tmp;
   double abskey;
   int j, itmp, first, mid, last;
   first = 0;
   last = n-1;
   if (ncut<first || ncut>last) return 0;
/* outer loop -- while mid != ncut */
label1:
   mid = first;
   abskey = cabs(a[mid]);
  for (j=first+1; j<=last; j++) {
     if (cabs(a[j]) > abskey) {
	 tmp = a[++mid];
	 itmp = ind[mid];
	 a[mid] = a[j];
	 ind[mid] = ind[j];
	 a[j]  = tmp;
	 ind[j] = itmp;
      }
   }
/*-------------------- interchange */
   tmp = a[mid];
   a[mid] = a[first];
   a[first]  = tmp;
   itmp = ind[mid];
   ind[mid] = ind[first];
   ind[first] = itmp;
/*-------------------- test for while loop */
   if (mid == ncut) return 0;
   if (mid > ncut) 
      last = mid-1;
   else
      first = mid+1;
   goto label1;
}
/*--------------- end of zqsplitC ----------------------------------------


|---------------------------------------------------------------------*/
int zSparTran(csptr amat, csptr bmat, int job, int flag)
{
/*----------------------------------------------------------------------
| Finds the transpose of a matrix stored in SpaFmt format.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SpaFmt format.
|
| job    = integer to indicate whether to fill the values (job.eq.1)
|          of the matrix (bmat) or only the pattern.
|
| flag   = integer to indicate whether the matrix has been filled
|          0 - no filled
|          1 - filled
|
| on return:
| ----------
| (bmat) = the transpose of (mata) stored in SpaFmt format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
  int i, j, *ind, pos, size=amat->n, *aja;
  complex double *ama=NULL;
  ind = (int *) Malloc(size*sizeof(int), "SparTran:1" );
  for (i=0; i<size; i++)
    ind[i] = 0;
  if(!flag) {
/*--------------------  compute lengths  */
    for (i=0; i<size; i++) {
      aja = amat->ja[i];
      for (j=0; j<amat->nzcount[i]; j++)
	ind[aja[j]]++;
    }
/*--------------------  allocate space  */
    for (i=0; i<size; i++) {
      bmat->ja[i] = (int *) Malloc(ind[i]*sizeof(int), "SparTran:2" );
      bmat->nzcount[i] = ind[i];
      if (job == 1) {
	bmat->ma[i] = (complex double *) Malloc(ind[i]*sizeof(complex double), "SparTran:3" );
      }
      ind[i] = 0;
    }
  }
/*--------------------  now do the actual copying  */
  for (i=0; i<size; i++) {
    aja = amat->ja[i];
    if (job == 1)
      ama = amat->ma[i];
    for (j=0; j<amat->nzcount[i]; j++) {
      pos = aja[j];
      bmat->ja[pos][ind[pos]] = i;
      if (job == 1)
	bmat->ma[pos][ind[pos]] = ama[j];
      ind[pos]++;
    }
  }
  free(ind);
  return 0;
}
/*-------------- end of SparTran ---------------------------------------
|---------------------------------------------------------------------*/

void zswapj(int v[], int i, int j)
{
   int temp;
   
   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
}
void zswapm(complex double v[], int i, int j)
{
   complex double temp;
   
   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
}
void swapm(double v[], int i, int j)
{
   double temp;
   
   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
}
  
void zqsortC(int *ja, complex double *ma, int left, int right, int abval)
/*----------------------------------------------------------------------
|
| qqsort: sort ma[left]...ma[right] into decreasing order
| from Kernighan & Ritchie
|
| ja holds the column indices
| abval = 1: consider absolute values
|         0: values
|
|---------------------------------------------------------------------*/
{
   int i, last;
   void zswapj(int *, int, int);
   void zswapm(complex double *, int, int);
   if (left >= right)  return;
   if (abval) {
      zswapj(ja, left, (left+right)/2);
      zswapm(ma, left, (left+right)/2);
      last = left;
      for (i=left+1; i<=right; i++) {
	 if (cabs(ma[i]) > cabs(ma[left])) {
	    zswapj(ja, ++last, i);
	    zswapm(ma, last, i);
	 }
      }
      zswapj(ja, left, last);
      zswapm(ma, left, last);
      zqsortC(ja, ma, left, last-1, abval);
      zqsortC(ja, ma, last+1, right, abval);
   }
   else {
      zswapj(ja, left, (left+right)/2);
      zswapm(ma, left, (left+right)/2);
      last = left;
      for (i=left+1; i<=right; i++) {
	 if (cabs(ma[i]) > cabs(ma[left])) {
	    zswapj(ja, ++last, i);
	    zswapm(ma, last, i);
	 }
      }
      zswapj(ja, left, last);
      zswapm(ma, left, last);
      zqsortC(ja, ma, left, last-1, abval);
      zqsortC(ja, ma, last+1, right, abval);
   }
}

int zroscalC(csptr mata, double *diag, int nrm)
{
/*---------------------------------------------------------------------
|
| This routine scales each row of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> row j is a zero row
|--------------------------------------------------------------------*/
/*--------------------   local variables    */
  int i, k;
  complex double *kr;
  double scal;
  for (i=0; i<mata->n; i++) {
    scal = 0.0;
    kr = mata->ma[i];
    if (nrm == 0) {
      for (k=0; k<mata->nzcount[i]; k++)
	if (cabs(kr[k]) > fabs(scal)) scal = cabs(kr[k]);
    }
    else if (nrm == 1) {
      for (k=0; k<mata->nzcount[i]; k++)
	scal += cabs(kr[k]);
    }
    else {  /* nrm = 2 */
      for (k=0; k<mata->nzcount[i]; k++)
	scal += cabs(kr[k]*kr[k]);
    }
    if (nrm == 2) scal = sqrt(scal);
    if (scal == 0.0) {
      scal = 1.0; 
      /* YS. return i+1; */
    }
    else 
      scal = 1.0 / scal;
    diag[i] = scal;
    for (k=0; k<mata->nzcount[i]; k++)
      kr[k] = kr[k] * scal;
  }
  return 0;
}
/*---------------end of roscalC-----------------------------------------
----------------------------------------------------------------------*/
int zcoscalC(csptr mata, double *diag, int nrm)
{
/*---------------------------------------------------------------------
|
| This routine scales each column of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SpaFmt form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> column j is a zero column
|--------------------------------------------------------------------*/
/*   local variables    */
   int i, j, k;
   complex double *kr;
   int *ki;
   for (i=0; i<mata->n; i++)
      diag[i] = 0.0;
/*---------------------------------------
|   compute the norm of each column
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      kr = mata->ma[i];
      ki = mata->ja[i];
      if (nrm == 0) {
	 for (k=0; k<mata->nzcount[i]; k++) {
	    j = ki[k];
	    if (cabs(kr[k]) > diag[j]) diag[j] = cabs(kr[k]);
	 }
      }
      else if (nrm == 1) {
         for (k=0; k<mata->nzcount[i]; k++)
            diag[ki[k]] += cabs(kr[k]);
      }
      else {  /*  nrm = 2 */
         for (k=0; k<mata->nzcount[i]; k++)
            diag[ki[k]] += cabs(kr[k]*kr[k]);
      }
   }
   if (nrm == 2) {
      for (i=0; i<mata->n; i++)
	 diag[i] = sqrt(diag[i]);
   }
/*---------------------------------------
|   invert
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      if (diag[i] == 0.0)
	/* return i+1;*/
	diag[i] = 1.0; 
      else 
	 diag[i] = 1.0 / diag[i];
   }
/*---------------------------------------
|   C = A * D
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      kr = mata->ma[i];
      ki = mata->ja[i];
      for (k=0; k<mata->nzcount[i]; k++)
	 kr[k] = kr[k] * diag[ki[k]];
   }
   return 0;
}
/*---------------end of coscalC-----------------------------------------
----------------------------------------------------------------------*/
void zdscale(int n, double *dd, complex double *x, complex double * y)
{ 
/* Computes  y == DD * x                               */
/* scales the vector x by the diagonal dd - output in y */
  int k;
  for (k=0; k<n; k++) 
    y[k] = dd[k]*x[k];
}

void zprintmat(FILE *ft, csptr A, int i0, int i1)
{
/*-------------------------------------------------------------+
| to dump rows i0 to i1 of matrix for debugging purposes       |
|--------------------------------------------------------------*/
  int i, k, nzi;
  int *row;
  complex double *rowm;
  for (i=i0; i<i1; i++)
    {
      nzi = A->nzcount[i];
      row = A->ja[i];
      rowm = A->ma[i];
      for (k=0; k< nzi; k++){
	fprintf(ft," row %d  a_real  %e a_imag %e ja %d \n", i+1, creal(rowm[k]), cimag(rowm[k]), row[k]+1);
	  }
    }
}

  
void qsortR2I(double *wa, int *cor1, int *cor2, int left, int right)
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into decreasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
{
   int i, last;
   void zswapj(int *, int, int);
   void swapm(double *, int, int);

   if (left >= right)  return;

   swapm(wa, left, (left+right)/2);
   zswapj(cor1, left, (left+right)/2);
   zswapj(cor2, left, (left+right)/2);
   last = left;
   for (i=left+1; i<=right; i++) {
      if (wa[i] > wa[left]) {
	 swapm(wa, ++last, i);
	 zswapj(cor1, last, i);
	 zswapj(cor2, last, i);
      }
   }
   swapm(wa, left, last);
   zswapj(cor1, left, last);
   zswapj(cor2, left, last);
   qsortR2I(wa, cor1, cor2, left, last-1);
   qsortR2I(wa, cor1, cor2, last+1, right);
}
void zqsort2C(int *ja, complex double *ma, int left, int right, int abval)
/*----------------------------------------------------------------------
|
| qqsort: sort ma[left]...ma[right] into increasing order
| from Kernighan & Ritchie
|
| ja holds the column indices
| abval = 1: consider absolute values
|         0: values
|
|---------------------------------------------------------------------*/
{
   int i, last;

   if (left >= right)  return;

   if (abval) {
      zswapj(ja, left, (left+right)/2);
      zswapm(ma, left, (left+right)/2);
      last = left;
      for (i=left+1; i<=right; i++) {
	 if (cabs(ma[i]) < cabs(ma[left])) {
	    zswapj(ja, ++last, i);
	    zswapm(ma, last, i);
	 }
      }
      zswapj(ja, left, last);
      zswapm(ma, left, last);
      zqsort2C(ja, ma, left, last-1, abval);
      zqsort2C(ja, ma, last+1, right, abval);
   }

   else {
      zswapj(ja, left, (left+right)/2);
      zswapm(ma, left, (left+right)/2);
      last = left;
      for (i=left+1; i<=right; i++) {
	 if (cabs(ma[i]) < cabs(ma[left])) {
	    zswapj(ja, ++last, i);
	    zswapm(ma, last, i);
	 }
      }
      zswapj(ja, left, last);
      zswapm(ma, left, last);
      zqsort2C(ja, ma, left, last-1, abval);
      zqsort2C(ja, ma, last+1, right, abval);
   }
}



void zqqsort(int *ja, complex double *ma, int left, int right)
/*----------------------------------------------------------------------
|
| qqsort: sort ja[left]...ja[right] into increasing order
| from Kernighan & Ritchie
|
| ma holds the real values
|
|---------------------------------------------------------------------*/
{
   int i, last;

   if (left >= right)  return;

   zswapj(ja, left, (left+right)/2);
   zswapm(ma, left, (left+right)/2);
   last = left;
   for (i=left+1; i<=right; i++) {
      if (ja[i] < ja[left]) {
	 zswapj(ja, ++last, i);
	 zswapm(ma, last, i);
      }
   }
   zswapj(ja, left, last);
   zswapm(ma, left, last);
   zqqsort(ja, ma, left, last-1);
   zqqsort(ja, ma, last+1, right);
}

void zhilosort(csptr mat, int abval, int hilo)
{
/*----------------------------------------------------------------------
|
| This routine sorts the entries in each row of a matrix from hi to low.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (mat) = a matrix stored in SpaFmt format.
|
| abval =   1: use absolute values of entries
|           0: use values
|
| hilo  =   1: sort in decreasing order
|           0: sort in increasing order
|
|
| on return:
| ----------
| (mat) = (mat) where each row is sorted.
|
|---------------------------------------------------------------------*/

   int j, n=mat->n, *nnz=mat->nzcount;
   void zqsortC(int *, complex double *, int, int, int);
   void zqsort2C(int *, complex double *, int, int, int);

   if (hilo)
      for (j=0; j<n; j++)
	 zqsortC(mat->ja[j], mat->ma[j], 0, nnz[j]-1, abval);

   else
      for (j=0; j<n; j++)
	 zqsort2C(mat->ja[j], mat->ma[j], 0, nnz[j]-1, abval);

   return;
}
/*------- end of hilosort ----------------------------------------------
|---------------------------------------------------------------------*/


void zqsort3i(int *wa, int *cor1, int *cor2, int left, int right)
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into increasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
{
   int i, last;
   void zswapj(int *, int, int);

   if (left >= right)  return;

   zswapj(wa, left, (left+right)/2);
   zswapj(cor1, left, (left+right)/2);
   zswapj(cor2, left, (left+right)/2);
   last = left;
   for (i=left+1; i<=right; i++) {
      if (wa[i] < wa[left]) {
	 zswapj(wa, ++last, i);
	 zswapj(cor1, last, i);
	 zswapj(cor2, last, i);
      }
   }
   zswapj(wa, left, last);
   zswapj(cor1, left, last);
   zswapj(cor2, left, last);
   zqsort3i(wa, cor1, cor2, left, last-1);
   zqsort3i(wa, cor1, cor2, last+1, right);
}

