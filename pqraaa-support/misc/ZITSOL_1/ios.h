#define MAX_LINE        256
#define MAX_MAT	   100
#define MAX_HBNAME      64

typedef struct _io_t {
  FILE *fout;                 /* output file handle              */
  char outfile[MAX_LINE];     /* output filename                 */
  char Fname[MAX_LINE];       /* matrix filename                 */
  char HBnameF[MAX_HBNAME];   /* HB name                         */
  char PrecMeth[MAX_LINE];    /* preconditioner being tested     */
  char type[4];               /* HB type                         */
  int ndim;                   /* matrix size                     */
  int nnz;                    /* number of nonzero               */
/* parameters from inputs -----------------------------------------*/
  int im;                     /* Dim of Krylov subspace [fgmr]   */
  int maxits;                 /* maximum number of fgmres iters  */
  double tol;                 /* tolerance for stopping fgmres   */
  int nparam;         /* number of tests for each preconditioner */
  int lfil0;                  /* initial lfil                    */
  int lfilInc;                /* increment for lfil              */
  double tol0;                /* initial drop tolerance          */
  double tolMul;              /* multiplier for tol              */    
  int fill_lev;               /* initial level of fill for ILUK  */
  int fill_lev_inc;               /* initial level of fill for ILUK  */
  int perm_type;              /* indset perms (0) or PQ perms (1) */
  int Bsize;                  /* block size - dual role. see input file */
                              /* for explanations */
/* deleted -- always take equal to one */
/* fill_lev_inc;           increment for level of fill for ILUK         */
/* result for output ----------------------------------------------*/
  double tm_p;                /* time for preconditioner (s)     */
  double tm_i;                /* time for iteration (s)          */
  double fillfact;                 /* memory used during precondition */
  int its;                    /* number of iterations            */
  double enorm;               /* error norm:          || x- x0|| */
  double rnorm;               /* final residual norm: ||Ax-Ax0|| */
} io_t;


/*-------------------- protos */
int zreadhb_c(int*, complex double**, int**, int**, io_t *pio, complex double **rhs,
	     complex double **guess, int*);
int zread_inputs( char *in_file, io_t *pio );
int zget_matrix_info( FILE *fmat, io_t *pio );

void zoutput_perm( int n, int *perm, FILE *f );
double sys_timer();
/*-------------------- end protos */
