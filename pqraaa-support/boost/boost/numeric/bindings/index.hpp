//
// Copyright (c) 2009 by Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_INDEX_HPP
#define BOOST_NUMERIC_BINDINGS_INDEX_HPP

#include <boost/numeric/bindings/rank.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <type_traits>
#include <cmath>

namespace boost {
namespace numeric {
namespace bindings {

template< typename T >
struct index_minor:
    std::conditional<
        is_column_major< T >::value,
        tag::index<1>,
        tag::index<
            std::max( tag::matrix::value, rank< T >::value )
        >
    >::type {};

template< typename T >
struct index_major:
    std::conditional<
        is_column_major< T >::value,
        tag::index<
            std::max( tag::matrix::value, rank< T >::value )
        >,
        tag::index<1>
    >::type {};


template< typename Index, typename TransTag >
struct index_trans {
    typedef Index type;
};

template<>
struct index_trans< tag::index<1>, tag::transpose > {
    typedef tag::index<2> type;
};

template<>
struct index_trans< tag::index<1>, tag::conjugate > {
    typedef tag::index<2> type;
};

template<>
struct index_trans< tag::index<2>, tag::transpose > {
    typedef tag::index<1> type;
};

template<>
struct index_trans< tag::index<2>, tag::conjugate > {
    typedef tag::index<1> type;
};


} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
