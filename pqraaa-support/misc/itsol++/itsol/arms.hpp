#ifndef ITSOL_ARMS_HPP
#define ITSOL_ARMS_HPP

#include <itsol/spar_mat.hpp>
#include <type_traits>
#include <vector>
#include <cassert>

/*---------------------------------------------------------------------
| MULTI-LEVEL BLOCK ILUT PRECONDITIONER.
| ealier  version:  June 23, 1999  BJS -- 
| version2: Dec. 07th, 2000, YS  [reorganized ]
| version 3 (latest) Aug. 2005.  [reorganized + includes ddpq]
+---------------------------------------------------------------------- 
| ON ENTRY:
| ========= 
| ( Amat ) = original matrix A stored in C-style Compressed Sparse
|            Row (CSR) format -- 
|            see LIB/heads.h for the formal definition of this format.
|
| ipar[0:17]  = integer array to store parameters for 
|       arms construction (arms2) 
|
|       ipar[0]:=nlev.  number of levels (reduction processes). 
|                       see also "on return" below. 
| 
|       ipar[1]:= level-reordering option to be used.  
|                 if ipar[1]==0 ARMS uses a block independent set ordering
|                  with a sort of reverse cutill Mc Kee ordering to build 
|                  the blocks. This yields a symmetric ordering. 
|                  in this case the reordering routine called is indsetC
|                 if ipar[1] == 1, then a nonsymmetric ordering is used.
|                  In this case, the B matrix is constructed to be as
|                  diagonally dominant as possible and as sparse as possble.
|                  in this case the reordering routine called is ddPQ.
|                 
|       ipar[2]:=bsize. for indset  Dimension of the blocks. 
|                  bsize is only a target block size. The size of 
|                  each block can vary and is >= bsize. 
|                  for ddPQ: this is only the smallest size of the 
|                  last level. arms will stop when either the number 
|                  of levels reaches nlev (ipar[0]) or the size of the
|                  next level (C block) is less than bsize.
|
|       ipar[3]:=iout   if (iout > 0) statistics on the run are 
|                       printed to FILE *ft
|
|       ipar[4-9] NOT used [reserved for later use] - set to zero.
| 
| The following set method options for arms2. Their default values can
| all be set to zero if desired. 
|
|       ipar[10-13] == meth[0:3] = method flags for interlevel blocks
|       ipar[14-17] == meth[0:3] = method flags for last level block - 
|       with the following meaning in both cases:
|            meth[0] nonsummetric permutations of  1: yes. affects rperm
|                    USED FOR LAST SCHUR COMPLEMENT 
|            meth[1] permutations of columns 0:no 1: yes. So far this is
|                    USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. 
|                    (so ipar[11] does no matter - enter zero). If 
|                    ipar[15] is one then ILUTP will be used instead 
|                    of ILUT. Permutation data stored in: perm2. 
|            meth[2] diag. row scaling. 0:no 1:yes. Data: D1
|            meth[3] diag. column scaling. 0:no 1:yes. Data: D2
|       all transformations related to parametres in meth[*] (permutation, 
|       scaling,..) are applied to the matrix before processing it 
| 
| ft       =  file for printing statistics on run
|
| droptol  = Threshold parameters for dropping elements in ILU 
|            factorization.
|            droptol[0:4] = related to the multilevel  block factorization
|            droptol[5:6] = related to ILU factorization of last block.
|            This flexibility is more than is really needed. one can use
|            a single parameter for all. it is preferable to use one value
|            for droptol[0:4] and another (smaller) for droptol[5:6]
|            droptol[0] = threshold for dropping  in L [B]. See piluNEW.c:
|            droptol[1] = threshold for dropping  in U [B].
|            droptol[2] = threshold for dropping  in L^{-1} F 
|            droptol[3] = threshold for dropping  in E U^{-1} 
|            droptol[4] = threshold for dropping  in Schur complement
|            droptol[5] = threshold for dropping  in L in last block
|              [see ilutpC.c]
|            droptol[6] = threshold for dropping  in U in last block
|              [see ilutpC.c]
|             This provides a rich selection - though in practice only 4
|             parameters are needed [which can be set to be the same 
              actually] -- indeed it makes sense to take
|             droptol[0] = droptol[1],  droptol[2] = droptol[3], 
|             and droptol[4] = droptol[5]
|
| lfil     = lfil[0:6] is an array containing the fill-in parameters.
|            similar explanations as above, namely:
|            lfil[0] = amount of fill-in kept  in L [B]. 
|            lfil[1] = amount of fill-in kept  in U [B].
|            lfil[2] = amount of fill-in kept  in E L\inv 
|            lfil[3] = amount of fill-in kept  in U \inv F
|            lfil[4] = amount of fill-in kept  in S    .
|            lfil[5] = amount of fill-in kept  in L_S  .
|            lfil[6] = amount of fill-in kept  in U_S 
|             
| tolind   = tolerance parameter used by the indset function. 
|            a row is not accepted into the independent set if 
|            the *relative* diagonal tolerance is below tolind.
|            see indset function for details. Good values are 
|            between 0.05 and 0.5 -- larger values tend to be better
|            for harder problems.
| 
| ON RETURN:
|=============
|
| (PreMat)  = arms data structure which consists of two parts:
|             levmat and ilsch. 
|
| ++(levmat)= permuted and sorted matrices for each level of the block 
|             factorization stored in PerMat4 struct. Specifically
|             each level in levmat contains the 4 matrices in:
|
|
|            |\         |       |
|            |  \   U   |       |
|            |    \     |   F   |
|            |  L   \   |       |
|            |        \ |       |
|            |----------+-------|
|            |          |       |
|            |    E     |       |
|            |          |       |
|            
|            plus a few other things. See LIB/heads.h for details.
|
| ++(ilsch) = IluSpar struct. If the block of the last level is:
|
|                        |  B    F |
|                  A_l = |         | 
|                        |  E    C |
|
|             then IluSpar contains the block C and an ILU
|             factorization (matrices L and U) for the last 
|             Schur complement [S ~ C - E inv(B) F ]
|             (modulo dropping) see LIB/heads.h for details.
|
| ipar[0]   = number of levels found (may differ from input value) 
|
*/

namespace itsol {

  template <typename T>
  class ARMS {
  } ;

  template <>
  class ARMS< double > {
    public:
      // nlevel:  number of levels (reduction processes). 
      // ordering: level-reordering option to be used.  
      //          if ordering==0 ARMS uses a block independent set ordering
      //           with a sort of reverse cutill Mc Kee ordering to build 
      //           the blocks. This yields a symmetric ordering. 
      //           in this case the reordering routine called is indsetC
      //          if ordering == 1, then a nonsymmetric ordering is used.
      //           In this case, the B matrix is constructed to be as
      //           diagonally dominant as possible and as sparse as possble.
      //           in this case the reordering routine called is ddPQ.
      //          
      // bsize. for indset  Dimension of the blocks. 
      //           bsize is only a target block size. The size of 
      //           each block can vary and is >= bsize. 
      //           for ddPQ: this is only the smallest size of the 
      //           last level. arms will stop when either the number 
      //           of levels reaches nlev (ipar[0]) or the size of the
      //           next level (C block) is less than bsize.
      //
      ARMS( int nlevel, int ordering, int bsize, double droptol, int level_fill, double tolind ) ;

      ~ARMS() ;

      template <typename M>
      int factorize( M& matrix ) {
        spar_mat_type< typename boost::numeric::bindings::value_type<M>::type > sp_mat ;
        spar_mat( matrix, sp_mat ) ;
        return this->factorize( sp_mat.n_, &sp_mat.nzcount_.front(), &sp_mat.ja_.front(), &sp_mat.ma_.front() ) ;
      }

      template <typename X>
      int solve( X& x ) {
        return solve( boost::numeric::bindings::size(x), boost::numeric::bindings::begin_value(x) ) ;
      }

    private:
      int factorize( int n, int *nzcount, int** ja, double** ma ) ;
      int solve( int n, double* x ) ;

    public:
      double droptol_[7] ;
      int    lfil_[7] ;
      int    ipar_[18] ;
      double tolind_ ;

    private:
      void* matrix_ ;
  } ;

} // namespace itsol

#endif
