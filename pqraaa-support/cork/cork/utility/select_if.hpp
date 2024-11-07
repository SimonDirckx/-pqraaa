#ifndef cork_utility_select_if_hpp
#define cork_utility_select_if_hpp

#include <type_traits>

namespace CORK {

  template <bool val>
  struct select_if {
  } ;

  template <>
  struct select_if< false > {
    template <typename A, typename B>
    static B const& apply( A const& a, B const& b ) {
      return b ;
    }
  } ;

  template <>
  struct select_if< true > {
    template <typename A, typename B>
    static A const& apply( A const& a, B const& b ) {
      return a ;
    }
  } ;

} // namespace CORK

#endif
