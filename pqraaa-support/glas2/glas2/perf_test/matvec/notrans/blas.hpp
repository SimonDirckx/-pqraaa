//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#define GLAS_BACKEND_BLAS
#define GLAS_BACKEND_DEFAULT backend_blas

#include "../../utils.hpp"
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>

template <typename T>
struct do_test {
void operator() ( int N, std::ostream& f ) const {
  N = int( sqrt( double(N) ) ) ;
  glas2::shared_vector< T > v( N ), w( N ) ;
  glas2::shared_matrix< T > a( N, N ) ;

  glas2::seed< decltype(std::abs(T())) > seed ;
  glas2::randomize( v, seed ); glas2::randomize( w, seed ) ;
  glas2::randomize( a, seed );

  int const n_loop( (32/sizeof(T)) * (100000000 / (N*N) ) ) ;
  std::cout << "Doing " << n_loop << " loops of square matvec() for size " << v.size() << std::endl ;
  boost::timer::cpu_timer timer ;
  for ( int i=0; i<n_loop; ++i ) {
    boost::numeric::bindings::blas::gemv( 1.0, a, v, 1.0, w ) ;
  }
  double time_elapsed = timer.elapsed().wall ;
  std::cout << "Sum : " << norm_2(w) << std::endl ;
  std::cout << "Elapsed time : " << time_elapsed << std::endl ;
  f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
}
} ;
