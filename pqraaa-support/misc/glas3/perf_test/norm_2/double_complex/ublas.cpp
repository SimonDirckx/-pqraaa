//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../ublas.hpp"

int main() {
  std::ofstream f( "ublas_time" ) ;
  utils::all_tests< do_test< std::complex<double> > >( f ) ;
  f.close() ;
  return 0 ;
}
