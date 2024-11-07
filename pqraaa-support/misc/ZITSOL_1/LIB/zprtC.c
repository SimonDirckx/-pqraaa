#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <complex.h>
#include "./zheads.h"
#ifdef _IBM
#define zprtmtc zprtmtc
#else
#define zprtmtc zprtmtc_ 
#endif 

/*---------------------------------------------------------------------------*
 * Complex version of prtC - Writes matrix in arms struct to HB file
 *---------------------------------------------------------------------------*/

int zprtC(csptr Amat, int io) 
{  
  int i, k, ko, nzi, n, nnz;
  int *ja = NULL, *ia = NULL, *col; 
  complex double  *aa = NULL, *rhs = NULL, *val;
  char guesol[3] = "NN";
  char title[73];
  char key[8]="ARMSMATx";
  char type[4] = "CUA";
  char fname[5]= "MAT"; 
  char ext[30] = "00";
  int ifmt = 6, job = 3;
/*-----------------------------------------------------------------------*/
  n   = Amat->n; 
  for (i=nnz=0; i<n; i++) 
    nnz += Amat->nzcount[i] ;
  ja  = (int *)malloc( sizeof(int) * nnz );
  aa  = (complex double *)malloc( sizeof(complex double) * nnz );
  ia  = (int *)malloc( sizeof(int) * (n+1));
  rhs = (complex double *)malloc( sizeof(complex double) * n );
/*-------------------- convert into a, ja, ia matrix                     */
  ko = 0;
  ia[0] = ko+1; 
  for (i=0; i<n; i++) {
    nzi = Amat->nzcount[i];
    val = Amat->ma[i] ;
    col = Amat->ja[i] ;
    for (k=0; k< nzi; k++) {
      aa[ko]  = val[k]; 
      ja[ko++]= col[k]+1; 
    }
    ia[i+1] = ko+1;
 
  }
/*-------------------- filename+string operations                        */
  /*  itoa(io,ext,10);*/
  sprintf(ext,"%02d",io);
  strncat(fname,ext,2); 
  strcpy(title,"Matrix from arms at level "); 
  strncat(title,ext,2); 
  sprintf(ext,", n =%5d",n);  
  strncat(title,ext,10); 
  for( i = strlen(title); i < 72; i++ ) title[i] = ' ';
  title[72] = 0;
  /* strcpy(key,"ARMSxx"); 
  for( i = strlen(key); i < 8; i++ ) key[i] = ' ';
  key[8] = 0;
  */
/*-------------------- do actual copying                               */
  zprtmtc( &n, &n, aa, ja, ia, rhs, guesol, title,
          key, type, &ifmt, &job, fname );
/*-------------------- done                                           */
  free(aa); 
  free (ia);
  free(ja);
  free(rhs);
  return 0;
}

