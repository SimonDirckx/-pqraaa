//
// Copyright (c) 2010 by Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_HAS_STATIC_SIZE_HPP
#define BOOST_NUMERIC_BINDINGS_HAS_STATIC_SIZE_HPP

//#include <boost/numeric/bindings/detail/property_map.hpp>
//#include <boost/numeric/bindings/tag.hpp>
//#include <boost/mpl/count_if.hpp>
//#include <boost/mpl/equal_to.hpp>
#include <utility>
#include <tuple>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T >
struct is_static_size_property: mpl::false_ {};

template< int N, int M >
struct is_static_size_property< std::pair< tag::size_type<N>, std::integral_constant<int,M> > >: std::true_type {};

} // namespace detail

template< typename T, typename Enable = void >
struct has_static_size: std::false_type {};

template< typename T >
struct has_static_size_one
: std::false_type
{} ;

template< typename T, typename ...Ts >
struct has_static_size_one< std::tuple< T, Ts... > >
: std::integral_constant< bool, has_static_size_one<std::tuple<Ts...>::value || detail::is_static_size_property< T >::value >
{} ;

template< typename T >
struct has_static_size<
        T,
        typename std::enable_if< detail::is_adaptable<T>::value >::type >:
has_static_size_one< typename detail::property_map_of< T >::type >
{} ;


/*
template< typename T >
struct has_static_size<
        T,
        typename std::enable_if< detail::is_adaptable<T>::value >::type >:

    // count the number of static size properties, 
    // should be equal to the rank of the object
    mpl::equal_to< 
        mpl::count_if<
            typename detail::property_map_of< T >::type,
            detail::is_static_size_property< mpl::_ >
        >,
        typename detail::property_at< T, tag::entity >::type
    >::type {};
*/
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
