//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_type_pass_argument_hpp
#define glas2_type_pass_argument_hpp

#include <type_traits>

namespace glas2 {

  template <typename Ref>
  class pass_argument {
    public:
      pass_argument( Ref ref )
      : ref_( ref )
      {}

      //typename std::add_reference<Ref>::type operator() { return ref_ ; }
      Ref operator()() { return ref_ ; }

    public:
      Ref ref_ ;
  } ;

} // namespace glas2

#endif
