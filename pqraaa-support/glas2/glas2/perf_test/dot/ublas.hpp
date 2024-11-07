//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_perf_test_dot_ublas_hpp
#define glas2_perf_test_dot_ublas_hpp

#include "../utils.hpp"
#include <glas2/bindings/vector.hpp>
#include <glas2/vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <fstream>

namespace ublas = boost::numeric::ublas ;

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
    glas2::randomize( glas2::bind_vector( v_ ), seed ); glas2::randomize( glas2::bind_vector( w_ ), seed ) ;
  }


  inline void run() {
    sum_ += ublas::inner_prod( v_, w_ ) ;
  }

  T                       sum_ ;
  ublas::vector< T >      v_ ;
  ublas::vector< T >      w_ ;
} ;

#endif
