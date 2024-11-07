//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#define GLAS_BACKEND_BLAS
#define GLAS_BACKEND_DEFAULT backend_blas

#include "../../utils.hpp"

#include <glas3/array/dense_array/algorithm/randomize.hpp>
#include <glas3/array/dense_array/algorithm/norm.hpp>
#include <glas3/bindings/array/dense_array/bind_vector.hpp>
#include <glas3/bindings/array/dense_array/bind_matrix.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/vector/adaptor.hpp>
#include <glas2/bindings/matrix/adaptor.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <string>

typedef typename boost::uniform_real<double>                        distribution_type ;
typedef typename boost::lagged_fibonacci607                         engine_type ;
typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
std::random_device rd ;
generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

#include <boost/timer.hpp>

template <typename T>
struct do_test {
void operator() ( int N, std::ostream& f ) const {
  N = int( sqrt(double(N)) ) ;
  glas2::shared_vector< T > v( N ), w( N ) ;
  glas2::shared_matrix< T, glas2::row_major > a( N, N ) ;

  glas3::randomize( glas3::bind_vector( v ), g ); glas3::randomize( glas3::bind_vector( w ), g ) ;
  glas3::randomize( glas3::bind_matrix( a ), g );

  int const n_loop( (32/sizeof(T)) * (100000000 / (N*N) ) ) ;
  std::cout << "Doing " << n_loop << " loops of square matvec() for size " << v.size() << std::endl ;
  boost::timer timer ;
  for ( int i=0; i<n_loop; ++i ) {
    boost::numeric::bindings::blas::gemv( 1.0, a, v, 1.0, w ) ;

  }
  double time_elapsed = timer.elapsed() ;
  std::cout << "Sum : " << glas3::norm_2( glas3::bind_vector( w ) ) << std::endl ;
  std::cout << "Elapsed time : " << time_elapsed << std::endl ;
  f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
}
} ;
