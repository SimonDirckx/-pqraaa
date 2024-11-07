//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_iterative_krylov_progress_report_hpp
#define glas_toolbox_iterative_krylov_progress_report_hpp

#include <boost/progress.hpp>
#include <string>

namespace glas {

  template <typename R>
  class progress_report {
    public:
      inline progress_report( std::ostream& s, R const& tol )
      : p_( int(-log(tol)), s )
      , tol_( tol )
      {}

    public:
      template <typename I>
      void operator() ( R const& r, I i ) {
        if (i==0) start_ = r ;
        double reference = -log(std::max(tol_, r/start_) ) - p_.count() ;
        if (reference>=0) p_ += reference ;
      }

    private:
      boost::progress_display p_ ;
      R                       tol_ ;
      R                       start_ ;
  } ;
}

#endif
