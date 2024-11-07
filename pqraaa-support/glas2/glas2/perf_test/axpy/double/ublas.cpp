//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../ublas.hpp"
#include <fstream>

int main() {
  std::ofstream f( "ublas_time" ) ;
  //utils::all_tests< do_test<double> >( f ) ;
  utils::all_loop_tests< do_one_loop<double> >( f ) ;
  f.close() ;
  return 0 ;
}
