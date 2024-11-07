//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_iterative_krylov_diagnostics_report_hpp
#define glas_toolbox_iterative_krylov_diagnostics_report_hpp

#include <boost/lexical_cast.hpp>
#include <string>

namespace glas {

  class diagnostics_report {
    public:
      inline diagnostics_report()
      : count_( 0 )
      {}

    public:
      template <typename R, typename I>
      inline void operator() ( R const& r, I i ) {
        count_ = i ;
        resid_ = r ;
      }

    public:
      inline std::string report() const {
        return "Number of iterations: " + boost::lexical_cast<std::string>(count_) + "\nResidual norm: " + boost::lexical_cast<std::string>(resid_) ;
      }

    private:
      int    count_ ;
      double resid_ ;
  } ;
}

#endif
