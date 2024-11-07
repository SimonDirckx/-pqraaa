#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../LIB/zheads.h"
#include "../LIB/zdefs.h" 
#include "../LIB/zprotos.h"  
#include "../ios.h" 
#include <sys/time.h>
#include <string.h>
#include <complex.h>

/*---------------------------------------------------------------------------*
 * main test driver for the zARMS2 preconditioner
 *---------------------------------------------------------------------------* 
 * This driver reads  input parameters from an input  file (default file
 * is "inputs"  - but one can  specify another file at  command line) --
 * some parameters are  also set `by hand' in  the function readin.  The
 * program  also  reads Harwell-Boeing  data  :  matrix  and right  hand
 * side(s).  If no right hand side is found, it constructs an artificial
 * one and finally solves the linear system using ARMS2. 
 *---------------------------------------------------------------------------*/

#define TOL_DD 0.8            /* diagonal dominance tolerance for ind-*/
                              /*-ependent sets                        */
/*-------------------- protos */
void output_header( io_t *pio );
void output_result( int lfil, io_t *pio, int iparam );
void set_arms_pars(io_t* io,  int Dscale, int *ipar, 
		   double *tolcoef, int *lfil);
int zreadhb_c(int *NN, complex double **AA, int **JA, int **IA, io_t *pio, 
	     complex double **rhs, complex double **guess, int *rsa);
int zread_inputs( char *in_file, io_t *pio );
int zget_matrix_info( FILE *fmat, io_t *pio );
void zrandvec (complex double *v, int n);
/*-------------------- end protos */

int main () {
  int ierr = 0;
/*-------------------------------------------------------------------
 * OPTIONS:
 * Plotting  = whether or not to dump gmres iteration infor (its, res)
 * diagscal  = diagonal scaling or not 
 *-----------------------------------------------------------------*/
  int plotting=0,  diagscal = 1;
  char pltfile[256];
  FILE *fits = NULL;

  double tol, tolind = TOL_DD;
  int j, lfil;   

  csptr csmat = NULL;    /* matrix in csr formt             */
  arms ArmsSt = NULL;    /* arms preconditioner structure   */
  SMatptr MAT = NULL;    /* Matrix structure for matvecs    */
  SPreptr PRE = NULL;    /* general precon  structure       */
  complex double *sol = NULL, *x = NULL, *rhs = NULL;
/*-------------------- method for incrementing lfil is set here */
  int lfil_arr[7]; 
  double droptol[7], dropcoef[7];
  int ipar[18];
/*-------------------- harwell boeing temp. arrays */
  complex double *AA;
  int *IA, *JA;
  int rsa; 
  int n; 
/*-------------------- IO-related */  
  FILE *flog = stdout;   /* to output stats */
  FILE *fmat = NULL;     /* matrix file     */
  io_t io;               /* structure for handling I/O functions +
                            other things */
  double tm1, tm2;
  int mat, numat, iparam, i;
  double terr;
  char line[MAX_LINE];
 
  MAT = (SMatptr)Malloc( sizeof(SMat), "main:MAT" );
  PRE = (SPreptr)Malloc( sizeof(SPre), "main:PRE" );
/*-------------------------------------------------------------------
 * reads matrix from Harwell Boeing file
 * solves using ARMS  preconditioned fgmres
 *-----------------------------------------------------------------*/
  memset(&io, 0, sizeof(io) );
/*-----------------------------------------------------------------*/
  if( zread_inputs( "inputs", &io ) != 0 ) {
    fprintf( flog, "ERROR reading inputs from file...\n" );
    goto ERROR_HANDLE;
  }
/*-----------------------------------------------------------------*/
  if( NULL == ( fmat = fopen( "matfile_hb", "r" ) ) ) {
    fprintf( flog, "Can't open matfile_hb...\n" );
    goto ERROR_HANDLE;
  }
  memset(line, 0, MAX_LINE );
  fgets(line, MAX_LINE, fmat );
  if( (numat = atoi( line ) ) <= 0 ) {
    fprintf( flog, "Invalid count of matrices...\n" );
    goto ERROR_HANDLE;
  }
/*-------------------- set parameters for arms */

  set_arms_pars(&io, diagscal, ipar, dropcoef, lfil_arr);

/*-------------------- open file OUT/ARMS.out for all performance
                       results of this run (all matrices and params) 
                       also set io->PrecMeth */
/* sprintf( io.outfile, "OUT/%s_ILUT.out", io.HBnameF );*/
    if (io.perm_type){
      strcpy(io.outfile,"OUT/ARMS_DDPQ.out");
      strcpy(io.PrecMeth,"ARMS_DDPQ");
    }
    else{
      strcpy(io.outfile,"OUT/ARMS.out");
      strcpy(io.PrecMeth,"ARMS");
    }
    if( NULL == ( io.fout = fopen( io.outfile, "w" ) ) ) {
      fprintf(flog,"Can't open output file %s...\n", io.outfile);
      goto ERROR_HANDLE;
    }
/*-------------------- LOOP through matrices -*/
  for( mat = 1; mat <= numat; mat++ ) {
    if( zget_matrix_info( fmat, &io ) != 0 ) {
      fprintf( flog, "Invalid format in matfile...\n" );
      goto ERROR_HANDLE;
    }
    fprintf( flog, "MATRIX: %s...\n", io.HBnameF );
/* Read in matrix and allocate memory------------------------------*/
    csmat = (csptr)Malloc( sizeof(zSparMat), "main:csmat" );
    ierr = zreadhb_c(&n,&AA, &JA, &IA, &io, &rhs, &sol, &rsa);
    if( ierr != 0 ) {
      fprintf( flog, "zreadhb_c error = %d\n", ierr );
      goto ERROR_HANDLE;
    }
/*-------------------- convert from CSR to SpaFmt matrix */
    if( ( ierr = zCSRcs( n, AA, JA, IA, csmat ) ) != 0 ) {
      fprintf( stderr, "ARMS: zCSRcs error\n" );
      return ierr;
    }
    free(IA); IA = NULL;
    free(AA); AA = NULL;
    free(JA); JA = NULL;
/*----------------------------------------------------------------------
|  The right-hand side is generated by assuming the solution is
|  a vector of ones.  To be changed if rhs is available from data.
|---------------------------------------------------------------------*/
    x   = (complex double *)Malloc( io.ndim * sizeof(complex double), "main" );
    for( i = 0; i < io.ndim; i++ ) 
      x[i] = 1.0 + 0.0*I;
    zmatvec(csmat, x, rhs) ;
    output_header(&io );
/*-------------------- set initial lfil and tol */ 
    lfil = io.lfil0;
    tol  = io.tol0;
/*-------------------- LOOP THROUGH PARAMETERS */
    for( iparam = 1; iparam <= io.nparam; iparam++ ) {
      fprintf( flog, "Parameter case  = %d\n", iparam );
      for (j=0; j<7; j++){
	lfil_arr[j] = lfil*((int) io.nnz/n); 
	droptol[j] = tol*dropcoef[j];
      }
      ArmsSt = (arms) Malloc(sizeof(zarmsMat),"main:ArmsSt");
      zsetup_arms(ArmsSt);
      fprintf( flog, "begin arms\n" );
      tm1 = sys_timer();
      ierr = zarms2(csmat, ipar, droptol, lfil_arr, tolind, ArmsSt, flog);
      tm2 = sys_timer();
      if( ierr != 0) { 
	fprintf( io.fout, " ** ARMS2 error - code %d...\n",ierr);
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      } 
      io.tm_p = tm2 - tm1;
      io.fillfact = (double)znnz_arms( ArmsSt, ipar[0], flog)/(double)(io.nnz + 1); 
      fprintf( flog, "ARMS ends, fill factor (mem used) = %f\n", io.fillfact );
      
/*-------------------- get rough idea of cond number */
      if(zcondestArms(ArmsSt, x, flog ) != 0 ) {
	fprintf( flog, "Not attempting iterative solution.\n" );
	fprintf( io.fout, "Not attempting iterative solution.\n" );
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }
/*-------------------- initial guess */
      zrandvec(x, n);
/* for(i=0; i < io.ndim; i++ ) 	x[i] = 0.0 + 0.0*I;        */
/*-------------------- create a file for printing
                       'its -- time -- res' info from fgmres */
     if(plotting ) { 
       sprintf( pltfile, "OUT/%s_ARMS_F%05d_T%08.6f", io.HBnameF, lfil,tol);
       if( NULL == ( fits = fopen( pltfile, "w" ) ) ) {
	 fprintf( flog, "Can't open output file %s...\n", pltfile );
	 goto ERROR_HANDLE;
       }
     } else 
       fits  =NULL;     
/*-------------------- set up the structs before calling fgmres */
      MAT->n = n;
      MAT->CSR = csmat;
      MAT->zmatvec = zmatvecCSR; 
      PRE->ARMS = ArmsSt;
      PRE->zprecon = zpreconARMS;
/*-------------------- call fgmr */
      io.its = io.maxits;
      tm1 = sys_timer();
      zfgmres(MAT, PRE, rhs, x, io.tol, io.im, &io.its, fits);
      tm2 = sys_timer();
      io.tm_i = tm2 - tm1;
      if( io.its < io.maxits ) 
 fprintf( flog, "param %03d OK: converged in %d steps...\n\n", iparam, io.its );
      else 
 fprintf( flog, "not converged in %d steps...\n\n", io.maxits );
       if( fits ) fclose( fits );
/*-------------------- calculate error norm */
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ ) {
	terr += pow(cabs( x[i] - 1.0 ), 2) ;
      }
      io.enorm = sqrt(terr);
      
/*-------------------- calculate residual norm */
      zmatvec(csmat, x, sol );
      terr = 0.0;
      for( i = 0; i < io.ndim; i++ )
	terr += pow(cabs( rhs[i] - sol[i] ),2);
      io.rnorm = sqrt(terr);
/*----------------------------- go to next param case */      
NEXT_PARA:
      output_result( lfil, &io, iparam );
      lfil += io.lfilInc;
      tol  *= io.tolMul;
      zcleanARMS( ArmsSt); 
    }
    
/*-------------------- NEXT_MAT: */
    zcleanCS(csmat);
    free( sol );
    free( x );
    free( rhs );
  }
  fclose( io.fout );
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  free(MAT);
  free(PRE);
  return 0;
  
ERROR_HANDLE:
  exit( -1 );
}
