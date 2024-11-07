//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_VALUE_OF_HPP
#define CORK_OPTIONS_VALUE_OF_HPP

#include <type_traits>
#include <tuple>

namespace CORK { namespace options {

  namespace detail {


    template <typename Key, typename Tuple, int Index, typename EnableIf=void>
    struct value_of_tuple
    {} ;


    template <typename Key, typename Tuple, int Index, typename EnableIf=void>
    struct value_of_tuple_index_larger_than_one
    {} ;

    template <typename Key, typename Tuple, int Index>
    struct value_of_tuple_index_larger_than_one< Key, Tuple, Index, typename std::enable_if< std::is_base_of<Key, typename std::tuple_element<Index-1,Tuple>::type>::value >::type >
    {
      static auto apply(Tuple const& tuple) {
        return std::get<Index-1>(tuple).value();
      }
      static auto& apply(Tuple& tuple) {
        return std::get<Index-1>(tuple).value() ;
      }
    } ;

    template <typename Key, typename Tuple, int Index>
    struct value_of_tuple_index_larger_than_one< Key, Tuple, Index, typename std::enable_if< !std::is_base_of<Key, typename std::tuple_element<Index-1,Tuple>::type>::value >::type >
    : value_of_tuple< Key, Tuple, Index-1>
    {} ;

    template <typename Key, typename Tuple, int Index>
    struct value_of_tuple< Key, Tuple, Index, typename std::enable_if< (Index>0) >::type >
    : value_of_tuple_index_larger_than_one< Key, Tuple, Index >
    {} ;

    template <typename Key, typename Tuple>
    struct value_of_tuple< Key, Tuple, 0 >
    {
      static auto apply(Tuple const& tuple) {
        return Key(tuple).value() ;
      }
      static auto& apply(Tuple& tuple) {
        return Key(tuple).value() ;
      }
    } ;

    template <typename Key, typename Tuple>
    struct value_of
    : value_of_tuple< Key, Tuple, std::tuple_size<Tuple>::value >
    {} ;

  } // namespace detail


  template <typename Key, typename Tuple>
  decltype(auto) value_of( Tuple const& tuple ) {
    return detail::value_of<Key,Tuple>::apply(tuple) ;
  }

  template <typename Key, typename Tuple>
  auto& value_of( Tuple& tuple ) {
    return detail::value_of<Key,Tuple>::apply(tuple) ;
  }

} } // namespace CORK::options

#endif
