//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_type_no_value_array_type
#define glas2_sparse_type_no_value_array_type

#include <glas2/sparse/type/no_value_type.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <type_traits>

namespace glas2 {

  struct no_value_type_array {
    typedef no_value_type value_type ;

    no_value_type_array() {}
    no_value_type_array( int n ) {}

    template <typename S>
    typename std::enable_if< std::is_integral<S>::value, no_value_type >::type operator()( S ) { return no_value_type() ; }

    template <typename S>
    typename vector_selection< no_value_type_array, S>::result_type operator()( S const& s ) { return vector_selection< no_value_type_array, S>::apply( no_value_type_array(), s ) ; }
  } ;

  template <typename S>
  struct vector_selection< no_value_type_array, S, typename std::enable_if< !std::is_integral<S>::value >::type > {
    typedef no_value_type_array result_type ;

    static result_type apply( no_value_type_array, S ) {
      return result_type() ;
    }
  } ;

} // namespace glas2


#endif
