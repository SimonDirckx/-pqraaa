//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_BLAS_OPTION_HPP
#define BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_BLAS_OPTION_HPP

#include <boost/numeric/bindings/tag.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace blas {
namespace detail {

template< typename Tag >
struct blas_option {};

template<>
struct blas_option< tag::transpose >: std::integral_constant< char, 'T' > {};

template<>
struct blas_option< tag::no_transpose >: std::integral_constant< char, 'N' > {};

template<>
struct blas_option< tag::conjugate >: std::integral_constant< char, 'C' > {};

template<>
struct blas_option< tag::upper >: std::integral_constant< char, 'U' > {};

template<>
struct blas_option< tag::lower >: std::integral_constant< char, 'L' > {};

template<>
struct blas_option< tag::unit >: std::integral_constant< char, 'U' > {};

template<>
struct blas_option< tag::non_unit >: std::integral_constant< char, 'N' > {};

template<>
struct blas_option< tag::left >: std::integral_constant< char, 'L' > {};

template<>
struct blas_option< tag::right >: std::integral_constant< char, 'R' > {};

} // namespace detail
} // namespace blas
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
