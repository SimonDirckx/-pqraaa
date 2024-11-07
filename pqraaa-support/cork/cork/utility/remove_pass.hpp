#ifndef cork_utility_remove_pass_hpp
#define cork_utility_remove_pass_hpp

namespace CORK {

  template <typename T>
  struct remove_pass_functor {
    typedef T type ;

    static type apply( T const& t ) { return t ; }
  } ;

  template <typename T>
  decltype (auto) remove_pass( T const& t) {
    return remove_pass_functor<T>::apply( t ) ;
  }

} // namespace CORK

#endif
