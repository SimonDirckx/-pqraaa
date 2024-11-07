//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_negate_hpp
#define glas2_concept_negate_hpp

namespace glas2 {

  struct negate {
    template <typename X>
    decltype( -X() ) operator() ( X const& x) const { return - x ; }
  } ;

}

#endif
