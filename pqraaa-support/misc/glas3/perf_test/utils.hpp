//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

//#include <glas2/utility/wall_clock_timer.hpp>
#include <boost/timer.hpp>
#include <fstream>
#include <iostream>

namespace utils {

template <class OneLoop>
struct do_loop {
  void operator() ( int N, std::ostream& f ) const {
    OneLoop one_loop( N ) ;
  
    int const n_loop( (32/sizeof(typename OneLoop::value_type)) * (100000000 / N ) ) ;
    std::cout << "Doing " << n_loop << " loops of " << one_loop.name() << " for problem size " << N << std::endl ;
    //glas::wall_clock_timer timer ;
    boost::timer timer ;
    for ( int i=0; i<n_loop; ++i ) {
      one_loop.run() ;
    }
    double time_elapsed = timer.elapsed() ;
    std::cout << "Elapsed time : " << time_elapsed << std::endl ;
    f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
  }
} ; // struct do_loop

template <class F>
void all_tests( std::ostream& f ) {
   F() ( 100, f ) ;
   F() ( 10000, f ) ;
   F() ( 1000000, f ) ;
}

template <class OneLoop>
void all_loop_tests( std::ostream& f ) {
  all_tests< do_loop<OneLoop> >( f ) ;
}

} // namespace utils
