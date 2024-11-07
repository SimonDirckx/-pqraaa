#ifndef cork_utility_tuple_hpp
#define cork_utility_tuple_hpp

#include <tuple>
#include <type_traits>

namespace CORK {

  namespace detail {
    template <int N>
    struct pop_front {
      template <typename T>
      static auto apply( T const& t ) {
        return std::tuple_cat( detail::pop_front<N-1>::apply(t), std::tuple( std::get<N>(t) ) ) ;
      }
    } ;

    template <>
    struct pop_front<1> {
      template <typename T>
      static auto apply( T const& t ) {
        return std::tuple( std::get<1>(t) ) ;
      }
    } ;

    template <>
    struct pop_front<0> {
      template <typename T>
      static auto apply( T const& t ) {
        return std::tuple<>() ;
      }
    } ;
  }

  template <class... T>
  auto pop_front( std::tuple<T...> const& t ) {
    return detail::pop_front<std::tuple_size<std::tuple<T...>>::value-1>::apply( t ) ;
  }

} // namespace CORK

#endif
