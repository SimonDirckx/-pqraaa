#include <itsol/iluk.hpp>

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

  ILUK<double>::ILUK( int level_fill )
  : lfil_( level_fill )
  , matrix_( (void*) new ILUSpar() )
  {}

  ILUK<double>::~ILUK()
  {
    cleanILU( static_cast<iluptr>(matrix_) ) ;
  }

  int ILUK<double>::factorize( int n, int *nzcount, int** ja, double** ma ) {
    SparMat csmat ;
    csmat.n = n ;
    csmat.nzcount = nzcount ;
    csmat.ja = ja ;
    csmat.ma = ma ;
    int info = ilukC( lfil_, &csmat, static_cast<iluptr>(matrix_), stdout );
    return info ;
  }

  int ILUK<double>::solve( int n, double*  x ) {
    assert( n == static_cast<iluptr>(matrix_)->n ) ;
    return lusolC( x, x,  static_cast<iluptr>(matrix_) );
  }

} // namespace itsol
