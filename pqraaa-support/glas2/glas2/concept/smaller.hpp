//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_smaller_hpp
#define glas2_concept_smaller_hpp

namespace glas2 {

  struct smaller {
    template <typename X, typename Y>
    bool operator() ( X const& x, Y const& y) const { return x < y ; }
  } ;

}

#endif
