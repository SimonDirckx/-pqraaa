//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_perf_test_dot_ublas_hpp
#define glas3_perf_test_dot_ublas_hpp

#include "../utils.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <fstream>

typedef typename boost::uniform_real<double>                        distribution_type ;
typedef typename boost::lagged_fibonacci607                         engine_type ;
typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
std::random_device rd ;
generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

namespace ublas = boost::numeric::ublas ;

template <typename T>
struct do_one_loop {
  typedef T value_type ;

  std::string name() const {
    return "axpy" ;
  }

  ~do_one_loop() {
    std::cout << "Sum : " << w_(0) << std::endl ;
  }

  do_one_loop( int N )
  : v_( N )
  , w_( N )
  , three_( 0.3 )
  {
	  for ( int i = 0; i < N; ++i ) {
		  v_(i) = g() ;
		  w_(i) = g() ;
	  }
  }


  inline void run() {
    w_.plus_assign( three_ * v_ ) ;
  }

  ublas::vector< T >      v_ ;
  ublas::vector< T >      w_ ;
  T const                 three_ ;
} ;

#endif
