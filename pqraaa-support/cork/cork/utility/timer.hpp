//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_timings_hpp
#define cork_timings_hpp

#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <chrono>

namespace CORK {

  class timer {
    public:
      typedef double value_type ;

    public:
      inline timer()
      {}

      inline void tic() {
        tic_ = std::chrono::system_clock::now() ;
      }

      inline value_type toc() {
        auto toc = std::chrono::system_clock::now() ;
        return std::chrono::duration<value_type >(toc - tic_).count() ;
      }

    private:
       decltype(std::chrono::system_clock::now()) tic_ ;
  } ;

} // namespace CORK


#endif
