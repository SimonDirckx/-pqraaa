//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../../utils.hpp"

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>
#include <glas3/array/dense_array/algorithm/randomize.hpp>
#include <glas3/array/dense_array/algorithm/multiply.hpp>
#include <glas3/array/dense_array/algorithm/norm.hpp>
#include <glas3/array/dense_array/algorithm/ops_assign.hpp>

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
  N = int( sqrt( double(N) ) ) ;
  glas3::dense_vector< T > v( glas3::no_init(), N ), w( glas3::no_init(), N ) ;
  glas3::dense_matrix< T > a( glas3::no_init(), {N, N} ) ;

  glas3::randomize( v, g ); glas3::randomize( w, g ) ;
  glas3::randomize( a, g );

  int const n_loop( (32/sizeof(T)) * (100000000 / (N*N) ) ) ;
  std::cout << "Doing " << n_loop << " loops of square matvec() for size " << v.size() << std::endl ;
  boost::timer timer ;
  for ( int i=0; i<n_loop; ++i ) {
    w += glas3::multiply( a, v ) ;
  }
  double time_elapsed = timer.elapsed() ;
  std::cout << "Sum : " << glas3::norm_2( w ) << std::endl ;
  std::cout << "Elapsed time : " << time_elapsed << std::endl ;
  f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
}
} ;
