//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_DETAIL_ADAPTOR_HPP
#define BOOST_NUMERIC_BINDINGS_DETAIL_ADAPTOR_HPP

#include <boost/numeric/bindings/is_numeric.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <type_traits>
#include <tuple>
#include <utility>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename Id, typename Enable = void >
struct adaptor {
/*    typedef mpl::map<
        mpl::pair< tag::value_type, void >
    > property_map;*/
   
};

template< typename T >
struct adaptor_access
: adaptor< typename std::remove_const<T>::type, T >
{};

/*
   template< typename T >
struct is_adaptable: is_numeric< typename mpl::at<
        typename adaptor< typename std::remove_const<T>::type, T >::property_map,
        tag::value_type >::type > {};

template< typename T, typename Enable = void >
struct adaptor_access {};

template< typename T >
struct adaptor_access< T, typename std::enable_if< is_adaptable<T>::value >::type >:
    adaptor< typename std::remove_const<T>::type, T > {};

*/

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#include <boost/numeric/bindings/detail/pod.hpp>
#include <boost/numeric/bindings/detail/property_map.hpp>

#endif
