#ifndef ITSOL_ILUK_HPP
#define ITSOL_ILUK_HPP

#include <itsol/spar_mat.hpp>
#include <type_traits>
#include <vector>
#include <cassert>

/*----------------------------------------------------------------------------
 * ILUK preconditioner
 * incomplete LU factorization with level of fill dropping
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * lofM     = level of fill: all entries with level of fill > lofM are
 *            dropped. Setting lofM = 0 gives BILU(0).
 * csmat    = matrix stored in SpaFmt format -- see heads.h for details
 *            on format
 * lu       = pointer to a ILUKSpar struct -- see heads.h for details
 *            on format
 * fp       = file pointer for error log ( might be stderr )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> error in lofC
 *            ierr  = -2  --> zero diagonal found
 * lu->n    = dimension of the matrix
 *   ->L    = L part -- stored in SpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in SpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonals of the input matrix must not be zero
 *--------------------------------------------------------------------------*/

namespace itsol {

  template <typename T>
  class ILUK {
  } ;

  template <>
  class ILUK< double > {
    public:
      ILUK( int level_fill ) ;

      ~ILUK() ;

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

    private:
      void* matrix_ ;
  } ;

} // namespace itsol

#endif
