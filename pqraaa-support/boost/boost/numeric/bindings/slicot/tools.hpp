//
// Copyright (c) 2011-
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_SLICOT_TOOLS_HPP
#define BOOST_NUMERIC_BINDINGS_SLICOT_TOOLS_HPP

#include <boost/numeric/bindings/is_column_major.hpp>

namespace boost { namespace numeric { namespace bindings { namespace slicot { namespace detail {

  template <typename Bool>
  struct trans_impl {
  } ;

  template <>
  struct trans_impl< boost::mpl::true_ > {
    BOOST_STATIC_CONSTANT( char, value='N') ;
  } ;

  template <>
  struct trans_impl< boost::mpl::false_ > {
    BOOST_STATIC_CONSTANT( char, value='T') ;
  } ;

  template <typename A>
  struct trans
  : trans_impl< typename is_column_major<A>::type >
  {} ;

} } } } } // namespace boost::numeric::bindings::slicot::detail

#endif
