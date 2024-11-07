//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_perf_test_dot_ublas_hpp
#define glas2_perf_test_dot_ublas_hpp

#include "../utils.hpp"
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <fstream>

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
    glas2::seed< decltype( std::abs(T() ) ) > seed ;
    glas2::randomize( glas2::bind_vector( v_ ), seed ); glas2::randomize( glas2::bind_vector( w_ ), seed ) ;
  }


  inline void run() {
    w_.plus_assign( three_ * v_ ) ;
  }

  ublas::vector< T >      v_ ;
  ublas::vector< T >      w_ ;
  T const                 three_ ;
} ;

#endif
