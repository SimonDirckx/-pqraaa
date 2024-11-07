//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_no_report_hpp
#define glas2_iterative_krylov_no_report_hpp

namespace glas2 { namespace iterative {

  class no_report {
    public:
      inline no_report()
      {}

    public:
      template <typename R, typename I>
      void operator() ( R const& , I ) {
      }
  } ;
} }

#endif
