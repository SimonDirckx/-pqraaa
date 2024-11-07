//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_DEFAULT_ORDER_HPP
#define BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_DEFAULT_ORDER_HPP

#include <boost/numeric/bindings/is_row_major.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <boost/numeric/bindings/detail/property_map.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace blas {
namespace detail {

template< typename T >
struct default_order {
    typedef typename std::conditional<
        bindings::detail::is_same_at< T, tag::value_transform, tag::conjugate >::value,
        typename std::conditional< is_row_major< T >::value, tag::column_major, tag::row_major >::type,
        tag::column_major
    >::type type;
};

} // namespace detail
} // namespace blas
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
