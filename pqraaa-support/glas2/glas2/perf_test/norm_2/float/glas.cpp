//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../glas.hpp"
#include <fstream>

int main() {
  std::ofstream f( "glas_time" ) ;
  utils::all_tests< do_test<float> >( f ) ;
  f.close() ;
  return 0 ;
}

