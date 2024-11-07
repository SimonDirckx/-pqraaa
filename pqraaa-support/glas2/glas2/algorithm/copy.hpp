#ifndef glas2_algorithm_copy_hpp
#define glas2_algorithm_copy_hpp

#include <glas2/type/copy_of.hpp>

namespace glas2 {

  template <typename E>
  copy_of<E> copy( E const& e ) {
    return copy_of<E>( e ) ;
  }

} // namespace glas2

#endif
