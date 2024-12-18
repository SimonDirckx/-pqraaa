//
// Copyright (c) 2009 by Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_HAS_LINEAR_ARRAY_HPP
#define BOOST_NUMERIC_BINDINGS_HAS_LINEAR_ARRAY_HPP

#include <boost/numeric/bindings/detail/property_map.hpp>
#include <boost/numeric/bindings/tag.hpp>

namespace boost {
namespace numeric {
namespace bindings {

template< typename T, typename Enable = void >
struct has_linear_array: std::false_type {};

template< typename T >
struct has_linear_array<
        T,
        typename std::enable_if< detail::is_adaptable<T>::value >::type >:
    detail::is_same_at< T, tag::data_structure, tag::linear_array > {};

} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
