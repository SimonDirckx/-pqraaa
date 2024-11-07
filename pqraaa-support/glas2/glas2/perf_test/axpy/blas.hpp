//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_perf_test_dot_blas_hpp
#define glas_perf_test_dot_blas_hpp

#include "../utils.hpp"
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/blas/level1.hpp>
#include <string>

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
    glas2::seed< decltype( std::abs(T() ) ) > seed ;
    glas2::randomize( v_, seed ); glas2::randomize( w_, seed ) ;
  }


  inline void run() {
    boost::numeric::bindings::blas::axpy( three_, v_, w_ ) ;
  }

  glas2::shared_vector< T > v_ ;
  glas2::shared_vector< T > w_ ;
  T const                   three_ ;
} ;

#endif
