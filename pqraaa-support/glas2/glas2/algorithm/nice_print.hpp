#ifndef glas2_algorithm_nice_print_hpp
#define glas2_algorithm_nice_print_hpp

#include <glas2/type/nice_print.hpp>

namespace glas2 {

  
  template <typename E>
  nice_print_type<E> nice_print( E const& e ) {return nice_print_type<E>( e ) ; }

} // namespace glas2

#endif
