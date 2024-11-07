//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_HAS_OPTION_HPP
#define CORK_OPTIONS_HAS_OPTION_HPP

#include <type_traits>
#include <tuple>

namespace CORK { namespace options {

  namespace detail {


    template <typename Key, typename Tuple, int Index, typename EnableIf=void>
    struct has_option_tuple
    {} ;


    template <typename Key, typename Tuple, int Index, typename EnableIf=void>
    struct has_option_tuple_index_larger_than_one
    {} ;

    template <typename Key, typename Tuple, int Index>
    struct has_option_tuple_index_larger_than_one< Key, Tuple, Index, typename std::enable_if< std::is_base_of<Key, typename std::tuple_element<Index-1,Tuple>::type>::value >::type >
    : std::true_type
    {} ;

    template <typename Key, typename Tuple, int Index>
    struct has_option_tuple_index_larger_than_one< Key, Tuple, Index, typename std::enable_if< !std::is_base_of<Key, typename std::tuple_element<Index-1,Tuple>::type>::value >::type >
    : has_option_tuple< Key, Tuple, Index-1>
    {} ;

    template <typename Key, typename Tuple, int Index>
    struct has_option_tuple< Key, Tuple, Index, typename std::enable_if< (Index>0) >::type >
    : has_option_tuple_index_larger_than_one< Key, Tuple, Index >
    {} ;

    template <typename Key, typename Tuple>
    struct has_option_tuple< Key, Tuple, 0 >
    : std::false_type
    {} ;

    template <typename Key, typename Tuple>
    struct has_option
    : has_option_tuple< Key, Tuple, std::tuple_size<Tuple>::value >
    {} ;

  } // namespace detail


  template <typename Key, typename Tuple>
  struct has_option
  : detail::has_option<Key,Tuple>
  {} ;

} } // namespace CORK::options

#endif
