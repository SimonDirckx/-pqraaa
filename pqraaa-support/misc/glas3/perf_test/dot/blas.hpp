//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_perf_test_dot_blas_hpp
#define glas_perf_test_dot_blas_hpp

#include "../utils.hpp"

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/algorithm/randomize.hpp>
#include <glas3/bindings/array/dense_array/adaptor_dense_vector.hpp>
#include <boost/numeric/bindings/blas/level1.hpp>

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

template <typename T>
struct do_one_loop {
  typedef T value_type ;

  std::string name() const {
    return "dot" ;
  }

  ~do_one_loop() {
    std::cout << "Sum : " << sum_ << std::endl ;
  }

  do_one_loop( int N )
  : sum_( 0.0 )
  , v_( glas3::no_init(), N )
  , w_( glas3::no_init(), N )
  {
    glas3::randomize( v_, g ); glas3::randomize( w_, g ) ;
  }


  inline void run() {
    sum_ += boost::numeric::bindings::blas::dot(v_, w_ ) ;
  }

  T                         sum_ ;
  glas3::dense_vector< T >  v_ ;
  glas3::dense_vector< T >  w_ ;
} ;

#endif
