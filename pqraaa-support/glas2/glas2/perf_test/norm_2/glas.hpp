//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_perf_test_norm_2_glas_hpp
#define glas_perf_test_norm_2_glas_hpp

#include "../utils.hpp"
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <fstream>

template <typename T>
struct do_test {
void operator() ( int N, std::ostream& f ) const {
  typedef T value_t ;
  glas2::shared_vector< value_t > v( N ) ;
  
  glas2::seed< decltype( std::abs(T() ) ) > seed ;
  glas2::randomize( v, seed );

  value_t sum( 0.0 ) ;
  int const n_loop( (32/sizeof(T)) * (100000000 / N ) ) ;
  std::cout << "Doing " << n_loop << " loops of norm_2() for vectors of length " << v.size() << std::endl ;
  boost::timer::cpu_timer timer ;
  for ( int i=0; i<n_loop; ++i ) {
    sum += value_t(norm_2(v)) ;
  }
  double time_elapsed = timer.elapsed().wall ;
  std::cout << "Sum : " << sum << std::endl ;
  std::cout << "Elapsed time : " << time_elapsed << std::endl ;
  f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
}
} ;

#endif
