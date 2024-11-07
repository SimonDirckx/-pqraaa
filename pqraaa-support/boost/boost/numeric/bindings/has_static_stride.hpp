//
// Copyright (c) 2010 by Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_HAS_STATIC_STRIDE_HPP
#define BOOST_NUMERIC_BINDINGS_HAS_STATIC_STRIDE_HPP

#include <boost/numeric/bindings/detail/property_map.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <utility>
#include <type_traits>
#include <tuple>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T >
struct is_static_stride_property: std::false_type {};

template< int N, int M >
struct is_static_stride_property< std::pair< tag::stride_type<N>, std::integral_constant<int, M> > >: std::true_type {};

} // namespace detail

template< typename T, typename Enable = void >
struct has_static_stride: std::false_type {};

template< typename T >
struct has_static_stride_count
: std::false_type
{} ;

template<typename T>
struct has_static_stride_count< std::tuple<T> >
: is_static_stride_property< T >
{} ;

template<typename T, typename ...Ts>
struct has_static_stride_count< std::tuple<T,Ts....> >
: std::integral_constant< is_static_stride_property< T >::value || has_static_stride_count< std::tuple<Ts....> >::value >
{} ;

template< typename T >
struct has_static_stride
        T,
        typename std::enable_if< detail::is_adaptable<T>::value >::type >
: has_static_stride_count< typename detail::property_map_of< T >::type >
{} ;

/*    // count the number of static stride properties, 
    // should be equal to the rank of the object
    mpl::equal_to< 
        mpl::count_if<
            typename detail::property_map_of< T >::type,
            detail::is_static_stride_property< mpl::_ >
        >,
        typename detail::property_at< T, tag::entity >::type
    >::type {};
*/
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
