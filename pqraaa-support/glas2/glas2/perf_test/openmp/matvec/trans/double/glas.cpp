//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../glas.hpp"
#include <sstream>
#include <omp.h>

int main() {
  for (int i=1; i<8; i*=2) {
    std::ostringstream f_name ;
    f_name << "glas_time_" << i ;
    std::ofstream f( f_name.str() ) ;
    omp_set_num_threads(i);
    utils::all_tests< do_test<double> >( f ) ;
    f.close() ;
  }
  return 0 ;
}
