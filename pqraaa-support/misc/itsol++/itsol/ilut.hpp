#ifndef ITSOL_ILUT_HPP
#define ITSOL_ILUT_HPP

#include <itsol/spar_mat.hpp>
#include <type_traits>
#include <vector>
#include <cassert>

/*----------------------------------------------------------------------------
 * ILUT preconditioner
 * incomplete LU factorization with dual truncation mechanism
 * NOTE : no pivoting implemented as yet in GE for diagonal elements
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * csmat    = matrix stored in SpaFmt format -- see heads.h for details
 * lu       = pointer to a ILUSpar struct -- see heads.h for details
 * lfil     = integer. The fill-in parameter. Each column of L and
 *            each column of U will have a maximum of lfil elements.
 *            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
 *            EARLIER VERSIONS. 
 *            lfil must be .ge. 0.
 * tol      = real*8. Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * fp       = file pointer for error log ( might be stdout )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> zero diagonal or zero col encountered
 * lu->n    = dimension of the matrix
 *   ->L    = L part -- stored in SpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in SpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonals of the input matrix must not be zero
 *----------------------------------------------------------------------------
 * Dual drop-off strategy works as follows. 
 *
 * 1) Theresholding in L and U as set by tol. Any element whose size
 *    is less than some tolerance (relative to the norm of current
 *    row in u) is dropped.
 *
 * 2) Keeping only the largest lfil elements in the i-th column of L
 *    and the largest lfil elements in the i-th column of U.
 *
 * Flexibility: one can use tol=0 to get a strategy based on keeping the
 * largest elements in each column of L and U. Taking tol .ne. 0 but lfil=n
 * will give the usual threshold strategy (however, fill-in is then
 * impredictible).
 *--------------------------------------------------------------------------*/

namespace itsol {

  template <typename T>
  class ILUT {
  } ;

  template <>
  class ILUT< double > {
    public:
      ILUT( int level_fill, double const& tol ) ;

      ~ILUT() ;

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
      int    lfil_ ;
      double tol_ ;

    private:
      void* matrix_ ;
  } ;

} // namespace itsol

#endif
