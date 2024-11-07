//  (C) Copyright Karl Meerbergen & Albert-Jan Yzelman 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_vector_concept_bsp_internal_vector_hpp
#define glas2_bsp_vector_concept_bsp_internal_vector_hpp

#include <glas2/vector/concept/vector.hpp>

namespace glas2 { namespace bsp {

  struct BSPInternalVector
  : glas2::Vector
  {
    typedef BSPInternalVector type ;
  } ;

  // Everything from Vector and, in addition,
  //
  // distribution_type distribution( )

} } // namespace glas::bsp

#endif
