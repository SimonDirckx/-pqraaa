#include <itsol/ilut.hpp>

extern "C" {
#include "globheads.h"
#include "defs.h" 
#include "protos.h"
#undef min
#undef max
}

#include <cassert>
#include <algorithm>

namespace itsol {

  ILUT<double>::ILUT( int level_fill, double const& tol )
  : lfil_( level_fill )
  , tol_( tol )
  , matrix_( (void*) new ILUSpar() )
  {}

  ILUT<double>::~ILUT()
  {
    cleanILU( static_cast<iluptr>(matrix_) ) ;
  }

  int ILUT<double>::factorize( int n, int *nzcount, int** ja, double** ma ) {
    SparMat csmat ;
    csmat.n = n ;
    csmat.nzcount = nzcount ;
    csmat.ja = ja ;
    csmat.ma = ma ;
    int info = ilut( &csmat, static_cast<iluptr>(matrix_), lfil_, tol_, stdout );
    return info ;
  }

  int ILUT<double>::solve( int n, double*  x ) {
    assert( n == static_cast<iluptr>(matrix_)->n ) ;
    return lusolC( x, x,  static_cast<iluptr>(matrix_) );
  }

} // namespace itsol
