#include <itsol/arms.hpp>

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

  ARMS<double>::ARMS( int nlevel, int ordering, int bsize, double droptol, int level_fill, double tolind )
  : tolind_( tolind )
  , matrix_( (void*) new armsMat() )
  {
    setup_arms( static_cast<arms>(matrix_) );

    std::fill( droptol_, droptol_+7, droptol ) ;
    std::fill( lfil_, lfil_+7, level_fill ) ;

    ipar_[0] = nlevel ;
    ipar_[1] = ordering ;
    ipar_[2] = bsize ;
    ipar_[3] = 0 ; // No printout

    ipar_[10] = ipar_[14] = 1 ; // rperm
    ipar_[11] = ipar_[15] = 1 ; // permutations of columns in ILUTP
    ipar_[12] = ipar_[16] = 1 ; // diagonal row scaling
    ipar_[13] = ipar_[17] = 1 ; // diagonal column scaling
  }

  ARMS<double>::~ARMS()
  {
    cleanARMS( static_cast<arms>(matrix_) ) ;
//    delete static_cast<arms>(matrix_) ;
  }

  int ARMS<double>::factorize( int n, int *nzcount, int** ja, double** ma ) {
    SparMat csmat ;
    csmat.n = n ;
    csmat.nzcount = nzcount ;
    csmat.ja = ja ;
    csmat.ma = ma ;
    int info = arms2( &csmat, ipar_, droptol_, lfil_, tolind_, static_cast<arms>(matrix_), stdout );
    return info ;
  }

  int ARMS<double>::solve( int n, double*  x ) {
    assert( n == static_cast<arms>(matrix_)->n ) ;
    return armsol2( x,  static_cast<arms>(matrix_) );
  }

} // namespace itsol
