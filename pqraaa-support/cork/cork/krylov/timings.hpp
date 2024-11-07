//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_timings_hpp
#define cork_krylov_timings_hpp

#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <chrono>

namespace CORK { namespace krylov {

  template <typename T>
  class timings {
    public:
      inline timings( int debug_level )
      : debug_level( debug_level )
      {
        reset() ;
        total_tic() ;
      }

      ~timings() {
        total_toc() ;
        if (debug_level>0) std::cout << *this << std::endl ;
      }

      void reset() {
        linearization_shift = 0.;
        linearization_solve = 0.;
        gram_schmidt = 0 ;
        other = 0 ;
        total = 0 ;
      }

      void total_tic() {
        total_tic_ = std::chrono::system_clock::now() ;
      }

      void total_toc() {
        auto toc = std::chrono::system_clock::now() ;
        total +=   std::chrono::duration<T>(toc - total_tic_).count() ;
      }

      void tic() {
        tic_ = std::chrono::system_clock::now() ;
      }

      void toc_linearization_shift() {
        auto toc = std::chrono::system_clock::now() ;
        linearization_shift += std::chrono::duration<T>(toc - tic_).count() ;
      }

      void toc_linearization_solve() {
        auto toc = std::chrono::system_clock::now() ;
        linearization_solve += std::chrono::duration<T>(toc - tic_).count() ;
      }

      void toc_gram_schmidt() {
        auto toc = std::chrono::system_clock::now() ;
        gram_schmidt += std::chrono::duration<T>(toc - tic_).count() ;
      }

      void toc_other() {
        auto toc = std::chrono::system_clock::now() ;
        other += std::chrono::duration<T>(toc - tic_).count() ;
      }

    private:
      typedef std::chrono::time_point<std::chrono::system_clock> time_stamp_type ;
      time_stamp_type total_tic_ ;
      time_stamp_type tic_ ;

    public:
      int debug_level ;
      T linearization_shift ;
      T linearization_solve ;
      T gram_schmidt ;
      T other ;
      T total ;
  } ;

  template <typename T>
  std::ostream& operator<<( std::ostream& os, timings<T> const& inf ) {
    os << "CORK timings:\n";
    os << "  new shift in linearization = " << inf.linearization_shift << "s.\n" ;
    os << "  solve with linearization = " << inf.linearization_solve << "s.\n" ;
    os << "  gram-schmidt orthogonalization = " << inf.gram_schmidt << "s.\n" ;
    os << "  other = " << inf.other << "s.\n" ;
    os << "  total = " << inf.total << "s.\n" ;
    return os ;
  }
   
} } // namespace CORK::krylov


#endif
