//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../blas.hpp"

int main() {
  std::ofstream f( "blas_time" ) ;
  utils::all_tests< do_test<double> >( f ) ;
  f.close() ;
  return 0 ;
}
