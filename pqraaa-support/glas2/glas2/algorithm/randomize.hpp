#ifndef glas2_algorithm_randomize_hpp
#define glas2_algorithm_randomize_hpp

#include <glas2/backend/current_backend.hpp>
#include <glas2/type/pass_reference.hpp>
#include <glas2/type/seed.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename S>
  auto randomize( V v, S& seed ) {
    return randomize( current_backend(), v, seed ) ;
  }

  template <typename V>
  auto randomize( V v, typename std::enable_if< pass_reference<V>::pass_by_value, int>::type i=0 ) {
    glas2::seed<double> seed ;
    return randomize( v, seed ) ;
  }

  template <typename V>
  auto randomize( V& v, typename std::enable_if< !pass_reference<V>::pass_by_value, int>::type i=0 ) {
    glas2::seed<double> seed ;
    return randomize( pass_reference<V>()(v), seed ) ;
  }

} // namespace glas2

#endif
