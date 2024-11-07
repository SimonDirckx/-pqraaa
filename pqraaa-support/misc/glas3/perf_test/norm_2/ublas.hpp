//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_perf_test_norm_2_ublas_hpp
#define glas_perf_test_norm_2_ublas_hpp

#include "../utils.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/timer.hpp>

typedef typename boost::uniform_real<double>                        distribution_type ;
typedef typename boost::lagged_fibonacci607                         engine_type ;
typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
std::random_device rd ;
generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

#include <boost/timer.hpp>
#include <fstream>

namespace ublas = boost::numeric::ublas ;

template <typename T>
struct do_test {
void operator() ( int N, std::ostream& f ) const {
  typedef T value_t ;
  ublas::vector< value_t > v( N ) ;
  
  for ( int i = 0; i < N; ++i ) {
	  v(i) = g() ;
  }

  value_t sum( 0.0 ) ;
  int const n_loop( (32/sizeof(T)) * (100000000 / N ) ) ;
  std::cout << "Doing " << n_loop << " loops of norm_2() for vectors of length " << v.size() << std::endl ;
  boost::timer timer ;
  for ( int i=0; i<n_loop; ++i ) {
    sum += norm_2(v) ;
  }
  double time_elapsed = timer.elapsed() ;
  std::cout << "Sum : " << sum << std::endl ;
  std::cout << "Elapsed time : " << time_elapsed << std::endl ;
  f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
}
} ;

#endif
