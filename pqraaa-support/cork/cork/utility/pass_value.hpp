#ifndef cork_utility_pass_value_hpp
#define cork_utility_pass_value_hpp

#include <cork/utility/remove_pass.hpp>

namespace CORK {

  template <typename T>
  struct pass_value_type {
    typedef T type ;

    pass_value_type( T const& t )
    : t_( t )
    {}

    T const& t_ ;
  } ;

  template <typename T>
  decltype (auto) pass_value( T const& t) {
    return pass_value_type<T>( t ) ;
  }

  template <typename T>
  struct remove_pass_functor< pass_value_type<T> >
  {
    typedef T type ;
    static T const& apply( pass_value_type<T> const& t ) { return t.t_ ; }
  } ;

} // namespace CORK

#endif
