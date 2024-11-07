//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_DETAIL_PROPERTY_MAP_HPP
#define BOOST_NUMERIC_BINDINGS_DETAIL_PROPERTY_MAP_HPP

#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/is_numeric.hpp>
#include <type_traits>
#include <tuple>
#include <utility>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename Tuple, typename Key >
struct property_map_has_key
: std::false_type
{} ;

template< typename T, typename... Ts, typename Key >
struct property_map_has_key< std::tuple< T, Ts... >, Key >
: std::conditional< std::is_same< Key, typename T::first_type >::value
                  , std::true_type
                  , property_map_has_key< std::tuple< Ts... >, Key >
                  >::type
                  {} ;


template< typename T, typename Key >
struct property_has_key: property_map_has_key< typename adaptor_access<T>::property_map, Key > {};

//template< typename T, typename Key >
//struct property_at {
//    typedef typename mpl::at< typename adaptor_access<T>::property_map, Key >::type type;
//};

template< typename T >
struct property_map_of {
    typedef typename adaptor_access<T>::property_map type;
};


template< typename T, typename Key >
struct property_map_at
: std::false_type
{} ;

template< typename T, typename... Ts, typename Key >
struct property_map_at< std::tuple< T, Ts... >, Key >
: std::conditional< std::is_same< Key, typename T::first_type >::value
                  , typename T::second_type
                  , typename property_map_at< std::tuple< Ts... >, Key >::type
                  >
                  {} ;

//template< typename Key >
//struct property_map_at< std::tuple<>, Key >
//: std::false_type
//{} ;


template< typename T, typename Key >
struct property_at
: property_map_at< typename property_map_of< T >::type, Key >
{
  //static_assert( property_has_key<T,Key>::value ) ;
};

template< typename T, typename Key, typename Value >
struct is_same_at: std::is_same< typename property_at< T, Key >::type, Value > {};


template< typename T >
struct is_adaptable: is_numeric< typename property_at< T, tag::value_type >::type >
{} ;


template< typename T, typename ...Ts >
struct property_insert {
  typedef typename property_map_of< T >::type properties;

  typedef decltype( std::tuple_cat( properties(), std::tuple<Ts...>() ) ) type ;
} ;

/*//
// Meta-function to insert multiple pairs into a map by using fold. Using the 
// provided mpl::insert can (only) insert elements in a map one-by-one.
//
template< typename T,
          typename P1 = mpl::na, typename P2 = mpl::na, typename P3 = mpl::na,
          typename P4 = mpl::na, typename P5 = mpl::na, typename P6 = mpl::na,
          typename P7 = mpl::na, typename P8 = mpl::na, typename P9 = mpl::na,
          typename P10 = mpl::na >
struct property_insert {

    typedef mpl::vector< P1, P2, P3, P4, P5, P6, P7, P8, P9, P10 > pair_vector;
    typedef typename property_map_of< T >::type properties;

    typedef typename mpl::fold<
        pair_vector,
        properties,
        mpl::if_<
            is_same< mpl::_2, mpl::void_ >,
            mpl::_1,
            mpl::insert< mpl::_1, mpl::_2 >
        >
    >::type type;

};*/

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
