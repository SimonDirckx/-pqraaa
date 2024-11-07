#ifndef cork_utility_ref_hpp
#define cork_utility_ref_hpp

#include <functional>
#include <memory>
#include <type_traits>

namespace CORK {

  namespace detail {

    template <typename T>
    struct deref {
      typedef T type ;
      static T& apply( T&  t ) { return t ;}
      static T const& apply( T const&  t ) { return t ;}
    } ;

    template <typename T>
    struct deref< std::reference_wrapper<T> > {
      typedef typename deref< typename std::remove_const<T>::type >::type type ;
      static auto& apply( std::reference_wrapper<T> const& t ) { return deref< typename std::remove_const<T>::type >::apply(t.get()) ;}
    } ;

    template <typename T>
    struct deref< T* > {
      typedef typename deref< typename std::remove_const<T>::type >::type type ;
      static auto& apply( T* t ) { return deref< typename std::remove_const<T>::type >::apply(*t) ;}
    } ;

    template <typename T>
    struct deref< std::unique_ptr<T> > {
      typedef typename deref< typename std::remove_const<T>::type >::type type ;
      static auto& apply( std::unique_ptr<T>& t ) { return deref< typename std::remove_const<T>::type >::apply(*t) ;}
      static auto& apply( std::unique_ptr<T> const& t ) { return deref< typename std::remove_const<T>::type >::apply(*t) ;}
    } ;

    template <typename T>
    struct deref< std::shared_ptr<T> > {
      typedef typename deref< typename std::remove_const<T>::type >::type type ;
      static auto& apply( std::shared_ptr<T>& t ) { return deref< typename std::remove_const<T>::type >::apply(*t) ;}
      static auto& apply( std::shared_ptr<T> const& t ) { return deref< typename std::remove_const<T>::type >::apply(*t) ;}
    } ;

  } // namespace detail

  template <typename T>
  auto& deref( T& t ) {
    return detail::deref<typename std::remove_const<T>::type>::apply(t) ;
  }

  template <typename T>
  struct deref_type
  : detail::deref<typename std::decay<T>::type>
  {} ;

} // namespace CORK

#endif
