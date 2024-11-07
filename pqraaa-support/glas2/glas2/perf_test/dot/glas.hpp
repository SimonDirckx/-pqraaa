//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_perf_test_dot_glas_hpp
#define glas2_perf_test_dot_glas_hpp

#include "../utils.hpp"
#include <glas2/vector.hpp>
#include <string>

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
  , v_( N )
  , w_( N )
  {
    glas2::seed< decltype( std::abs(T() ) ) > seed ;
    glas2::randomize( v_, seed ); glas2::randomize( w_, seed ) ;
  }


  inline void run() {
    sum_ += T(glas2::inner_prod( v_, w_ )) ;
  }

  T                       sum_ ;
  glas2::shared_vector< T > v_ ;
  glas2::shared_vector< T > w_ ;
} ;

#endif
