#ifndef cork_utility_pass_lvalue_reference_hpp
#define cork_utility_pass_lvalue_reference_hpp

#include <cork/utility/remove_pass.hpp>

namespace CORK {

  template <typename T>
  struct pass_lvalue_reference_type {
    typedef T const& type ;

    pass_lvalue_reference_type( T const& t )
    : t_( t )
    {}

    T const& t_ ;
  } ;

  template <typename T>
  decltype (auto) pass_lvalue_reference( T const& t) {
    return pass_lvalue_reference_type<T>( t ) ;
  }

  template <typename T>
  struct remove_pass_functor< pass_lvalue_reference_type<T> >
  {
    typedef T const& type ;
    static T const& apply( pass_lvalue_reference_type<T> const& t ) { return t.t_ ; }
  } ;

} // namespace CORK

#endif
