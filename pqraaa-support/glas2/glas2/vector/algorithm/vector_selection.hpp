#ifndef glas2_vector_algorithm_vector_selection_hpp
#define glas2_vector_algorithm_vector_selection_hpp

#include <glas2/vector/type/all.hpp>

namespace glas2 {

  template <typename V, typename S, typename EnableIf=void>
  struct vector_selection {
  } ;

/*  template <typename V>
  struct vector_selection< V, all > {
    typedef V result_type ;
   
    static result_type apply( V v, all ) {
      return v ;
    }
  } ;*/

}

#endif
