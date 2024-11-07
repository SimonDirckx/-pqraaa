//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_concept_need_wrapper_hpp
#define glas2_concept_need_wrapper_hpp

#include <boost/mpl/bool.hpp>

namespace glas2 {

  template <typename X, typename EnableIf=void>
  struct need_wrapper
  : boost::mpl::true_
  {} ;
}

#endif
