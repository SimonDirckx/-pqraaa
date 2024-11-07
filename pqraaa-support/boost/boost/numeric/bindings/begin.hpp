//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_BEGIN_HPP
#define BOOST_NUMERIC_BINDINGS_BEGIN_HPP

#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/linear_iterator.hpp>
#include <boost/numeric/bindings/rank.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <cassert>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {

namespace detail {

template< typename T, typename Tag >
struct begin_impl {};

template< typename T >
struct begin_impl< T, tag::value > {
    typedef typename bindings::value_type< T>::type* result_type;

    static result_type invoke( T& t ) {
        return adaptor_access<T>::begin_value( t );
    }
};

template< typename T, int Dimension >
struct begin_impl<T, tag::addressing_index<Dimension> > {

    typedef tag::addressing_index<Dimension> tag_type;

    typedef linear_iterator<
        typename bindings::value_type< T>::type,
        typename result_of::stride< T, tag_type >::type
    > result_type;

    static result_type invoke( T& t ) {
        return result_type( begin_value( t ), stride(t, tag_type() ) );
    }
};

template< typename T >
struct begin_impl< T, tag::index_major > {
    typedef typename property_at< T, tag::index_type >::type* result_type;

    static result_type invoke( T& t ) {
        return adaptor_access<T>::begin_index_major( t );
    }
};

template< typename T >
struct begin_impl< T, tag::compressed_index_major > {
    typedef typename property_at< T, tag::index_type >::type* result_type;

    static result_type invoke( T& t ) {
        return adaptor_access<T>::begin_compressed_index_major( t );
    }
};

template< typename T >
struct begin_impl< T, tag::index_minor > {
    typedef typename property_at< T, tag::index_type >::type* result_type;

    static result_type invoke( T& t ) {
        return adaptor_access<T>::begin_index_minor( t );
    }
};

} // namespace detail

namespace result_of {

template< typename T, typename Tag = tag::addressing_index<1> >
struct begin {
    static_assert( (is_tag<Tag>::value) );
    typedef typename detail::begin_impl<T,Tag>::result_type type;
};

} // namespace result_of

//
// Free Functions
//

//
// Overloads like begin( t, tag )
//
template< typename T, typename Tag >
inline typename result_of::begin<T,Tag>::type
begin( T& t, Tag ) {
    return detail::begin_impl<T,Tag>::invoke( t );
}

template< typename T, typename Tag >
inline typename result_of::begin<const T,Tag>::type
begin( const T& t, Tag ) {
    return detail::begin_impl<const T,Tag>::invoke( t );
}

// Overloads for types with rank <= 1 (scalars, vectors)
// In theory, we could provide overloads for matrices here, too, 
// if their minimal_rank is at most 1.

template< typename T >
//typename boost::enable_if< mpl::less< rank<T>, mpl::int_<2> >,
typename std::enable_if< rank<T>::value < 2,
    typename result_of::begin< T >::type >::type
begin( T& t ) {
    return detail::begin_impl< T, tag::addressing_index<1> >::invoke( t );
}

template< typename T >
//typename boost::enable_if< mpl::less< rank<T>, mpl::int_<2> >,
typename std::enable_if< rank<T>::value < 2,
    typename result_of::begin< const T >::type >::type
begin( const T& t ) {
    return detail::begin_impl< const T, tag::addressing_index<1> >::invoke( t );
}

namespace result_of {
  template< typename T >
  struct begin1 {
      typedef typename detail::begin_impl< T, tag::addressing_index<1> >::result_type type ;
  } ;
  template< typename T >
  struct begin2 {
      typedef typename detail::begin_impl< T, tag::addressing_index<2> >::result_type type ;
  } ;
  template< typename T >
  struct begin3 {
    typedef typename detail::begin_impl< T, tag::addressing_index<3> >::result_type type ;
  } ;
  template< typename T >
  struct begin_value {
    typedef typename detail::begin_impl< T, tag::value >::result_type type ;
  } ;
  template< typename T >
  struct begin_row {
    typedef typename detail::begin_impl< T, tag::addressing_index<1> >::result_type type ;
  } ;
  template< typename T >
  struct begin_column {
    typedef typename detail::begin_impl< T, tag::addressing_index<2> >::result_type type ;
  } ;
  template< typename T >
  struct begin_index_major {
    typedef typename detail::begin_impl< T, tag::index_major >::result_type type ;
  } ;
  template< typename T >
  struct begin_compressed_index_major {
    typedef typename detail::begin_impl< T, tag::compressed_index_major >::result_type type ;
  } ;
  template< typename T >
  struct begin_index_minor {
    typedef typename detail::begin_impl< T, tag::index_minor >::result_type type ;
  } ;
}

template< typename T >
typename result_of::begin1<T>::type begin1( T& t) {
  return detail::begin_impl< T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin2<T>::type begin2( T& t) {
  return detail::begin_impl< T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin3<T>::type begin3( T& t) {
  return detail::begin_impl< T, tag::addressing_index<3> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_value<T>::type begin_value( T& t) {
  return detail::begin_impl< T, tag::value >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_row<T>::type begin_row( T& t) {
  return detail::begin_impl< T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_column<T>::type begin_column( T& t) {
  return detail::begin_impl< T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_index_major<T>::type begin_index_major( T& t) {
  return detail::begin_impl< T, tag::index_major >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_compressed_index_major<T>::type begin_compressed_index_major( T& t) {
  return detail::begin_impl< T, tag::compressed_index_major >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_index_minor<T>::type begin_index_minor( T& t) {
  return detail::begin_impl< T, tag::index_minor >::invoke( t ) ;
}

template< typename T >
typename result_of::begin1<const T>::type begin1( const T& t) {
  return detail::begin_impl< const T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin2<const T>::type begin2( const T& t) {
  return detail::begin_impl< const T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin3<const T>::type begin3( const T& t) {
  return detail::begin_impl< const T, tag::addressing_index<3> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_value<const T>::type begin_value( const T& t) {
  return detail::begin_impl< const T, tag::value >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_row<const T>::type begin_row( const T& t) {
  return detail::begin_impl< const T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_column<const T>::type begin_column( const T& t) {
  return detail::begin_impl< const T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_index_major<const T>::type begin_index_major( const T& t) {
  return detail::begin_impl< const T, tag::index_major >::invoke( t ) ;
}
template< typename T >
typename result_of::begin<const T>::type begin_compressed_index_major( const T& t) {
  return detail::begin_impl< const T, tag::compressed_index_major >::invoke( t ) ;
}
template< typename T >
typename result_of::begin_index_minor<const T>::type begin_index_minor( const T& t) {
  return detail::begin_impl< const T, tag::index_minor >::invoke( t ) ;
}

/*
#define GENERATE_BEGIN_INDEX( z, which, unused ) \
GENERATE_FUNCTIONS( begin, which, tag::addressing_index<which> )

GENERATE_FUNCTIONS( begin, _value, tag::value )
GENERATE_FUNCTIONS( begin, _row, tag::addressing_index<1> )
GENERATE_FUNCTIONS( begin, _column, tag::addressing_index<2> )

GENERATE_FUNCTIONS( begin, _index_major, tag::index_major )
GENERATE_FUNCTIONS( begin, _compressed_index_major, tag::compressed_index_major )
GENERATE_FUNCTIONS( begin, _index_minor, tag::index_minor )
*/
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
