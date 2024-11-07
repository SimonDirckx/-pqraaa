//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_perf_test_dot_glas_hpp
#define glas2_perf_test_dot_glas_hpp

#include "../utils.hpp"
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <string>

template <typename T>
struct do_one_loop {
  typedef T value_type ;

  std::string name() const {
    return "axpy" ;
  }

  do_one_loop( int N )
  : v_( N )
  , w_( N )
  , vv_( v_ )
  , ww_( w_ )
  , three_( 0.3 )
  {
    glas2::seed< decltype( std::abs(T() ) ) > seed ;
    glas2::randomize( v_, seed ); glas2::randomize( w_, seed ) ;
  }

  ~do_one_loop() {
    std::cout << "Sum : " << w_(0) << std::endl ;
  }


  inline void run() {
    //auto e = glas2::binary_operation<T,glas2::contiguous_vector<T,std::ptrdiff_t>,glas2::plus >( three_, vv_) ;
//    w_ += glas2::binary_operation<T,glas2::contiguous_vector<T,std::ptrdiff_t>,glas2::plus >( three_, v_) ;
//    ww_ += three_ * vv_ ;
/*    for (int i=0; i<vv_.size(); ++i) {
      //ww_(i) += three_ * vv_(i) ;
      w_(i) += e(i) ;
    }*/
  //  glas2::plus_assign( glas2::current_backend(), ww_, three_*vv_ ) ;
    ww_ += three_ * vv_ ;
  }

  glas2::shared_vector< T > v_ ;
  glas2::shared_vector< T > w_ ;
  glas2::contiguous_vector< T, std::ptrdiff_t > vv_ ;
  glas2::contiguous_vector< T, std::ptrdiff_t > ww_ ;
  T const                  three_ ;
} ;

#endif
