//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_file_report_hpp
#define glas2_iterative_krylov_file_report_hpp

#include <iosfwd>
#include <string>

namespace glas2 { namespace iterative {

  class file_report {
    public:
      inline file_report( std::ostream& s )
      : s_( s )
      {}

    public:
      template <typename R, typename I>
      void operator() ( R const& r, I i ) {
        s_ << i << std::string(" ") << r << std::endl ;
      }

    private:
      std::ostream& s_ ;
  } ;
} }

#endif
