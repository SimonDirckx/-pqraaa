//
// Copyright (c) 2009 by Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_IS_NUMERIC_HPP
#define BOOST_NUMERIC_BINDINGS_IS_NUMERIC_HPP

#include <boost/numeric/bindings/is_real.hpp>
#include <boost/numeric/bindings/is_complex.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {

template< typename T >
struct is_numeric: std::integral_constant< bool, std::is_arithmetic<T>::value || is_complex<T>::value > {};

} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
