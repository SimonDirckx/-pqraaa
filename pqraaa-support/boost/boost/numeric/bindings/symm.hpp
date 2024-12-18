//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_SYMM_HPP
#define BOOST_NUMERIC_BINDINGS_SYMM_HPP

#include <boost/numeric/bindings/detail/basic_wrapper.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <utility>

namespace boost {
namespace numeric {
namespace bindings {
namespace result_of {

template< typename T >
struct symm {
    typedef detail::basic_wrapper<
        T,
        std::pair< tag::matrix_type, tag::symmetric >
    > type;
};

} // namespace result_of

template< typename T >
typename result_of::symm< T >::type const symm( T& underlying ) {
    return typename result_of::symm< T >::type( underlying );
}

template< typename T >
typename result_of::symm< const T >::type const symm( const T& underlying ) {
    return typename result_of::symm< const T >::type( underlying );
}

} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
