//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_concept_concept_hpp
#define glas3_concept_concept_hpp

namespace glas3 {

  template <typename X, typename EnableIf=void>
  struct concept
  {
    typedef void type ;
  } ;

}

#endif
