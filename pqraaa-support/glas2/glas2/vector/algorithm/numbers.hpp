#ifndef glas2_vector_algorithm_numbers_hpp
#define glas2_vector_algorithm_numbers_hpp

#include <glas2/vector/type/numbers_vector.hpp>
#include <initializer_list>

namespace glas2 {

  template <typename T>
  numbers_vector< std::initializer_list<T> > numbers( std::initializer_list<T> const& list ) {
    return numbers_vector< std::initializer_list<T> >( list ) ;
  }

} // namespace glas2

#endif
