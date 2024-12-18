//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_DETAIL_COPY_CONST_HPP
#define BOOST_NUMERIC_BINDINGS_DETAIL_COPY_CONST_HPP

#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename Source, typename Target >
struct copy_const
: std::conditional< std::is_const<Source>::value, typename  std::add_const<Target>::type, Target >
{} ;

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
