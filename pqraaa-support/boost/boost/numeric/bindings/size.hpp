//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_SIZE_HPP
#define BOOST_NUMERIC_BINDINGS_SIZE_HPP

#include <boost/numeric/bindings/detail/get.hpp>
#include <boost/numeric/bindings/rank.hpp>
#include <boost/numeric/bindings/addressing_index.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename AddressingIndex, typename Enable = void >
struct size_impl {

    typedef typename tag::size_type< AddressingIndex::value > key_type;
    typedef typename result_of_get< T, key_type >::type result_type;

    static result_type invoke( const T& t ) {
        return get< key_type >( t );
    }

};

template< typename T, typename AddressingIndex >
struct size_impl< T, AddressingIndex,
        typename std::enable_if< (AddressingIndex::value > rank<T>::value) && is_same_at< T, tag::size_type<1>, std::ptrdiff_t >::value
                               >::type > {
//        typename boost::enable_if< typename mpl::and_<
//            mpl::greater< AddressingIndex, rank<T> >,
//            is_same_at< T, tag::size_type<1>, std::ptrdiff_t >
//        >::type >::type > {

    typedef std::ptrdiff_t result_type;

    static result_type invoke( const T& t ) {
        return std::min< std::ptrdiff_t >( size_impl<T, tag::addressing_index<1> >::invoke(t), 1 );
    }

};

template< typename T, typename AddressingIndex >
struct size_impl< T, AddressingIndex,
        typename std::enable_if< (AddressingIndex::value > rank<T>::value) && !is_same_at< T, tag::size_type<1>, std::ptrdiff_t >::value
                               >::type > {
//        typename boost::enable_if< typename mpl::and_<
//            mpl::greater< AddressingIndex, rank<T> >,
//            mpl::not_< is_same_at< T, tag::size_type<1>, std::ptrdiff_t > >
//        >::type >::type > {

    typedef std::integral_constant< int, std::min( property_at< T, tag::size_type<1> >::value, 1 ) > result_type ;
//    typedef typename mpl::min<
//        typename property_at< T, tag::size_type<1> >::type,
//        mpl::int_<1>
//    >::type result_type;

    static result_type invoke( const T& t ) {
        return result_type();
    }

};

} // namespace detail


namespace result_of {

template< typename T, typename Tag = tag::addressing_index<1> >
struct size {
    static_assert( is_tag<Tag>::value );
    typedef typename detail::size_impl< T, Tag >::result_type type;
};

} // namespace result_of

//
// Overloads for free template functions size( x, tag ), 
//
template< typename T, typename Tag >
inline typename result_of::size< const T, Tag >::type
size( const T& t, Tag ) {
    return detail::size_impl< const T, Tag >::invoke( t );
}

// Overloads for free template function size( x )
// Valid for types with rank <= 1 (scalars, vectors)
// In theory, we could provide overloads for matrices here, too, 
// if their minimal_rank is at most 1.

template< typename T >
//typename boost::enable_if< mpl::less< rank<T>, mpl::int_<2> >,
typename std::enable_if< rank<T>::value < 2,
    typename result_of::size< const T >::type >::type
size( const T& t ) {
    return detail::size_impl< const T, tag::addressing_index<1> >::invoke( t );
}

namespace result_of {
  template< typename T >
  struct size1 {
      typedef typename detail::size_impl< T, tag::addressing_index<1> >::result_type type ;
  } ;
  template< typename T >
  struct size2 {
      typedef typename detail::size_impl< T, tag::addressing_index<2> >::result_type type ;
  } ;
  template< typename T >
  struct size3 {
      typedef typename detail::size_impl< T, tag::addressing_index<3> >::result_type type ;
  } ;
  template< typename T >
  struct size_row {
      typedef typename detail::size_impl< T, tag::addressing_index<1> >::result_type type ;
  } ;
  template< typename T >
  struct size_column {
      typedef typename detail::size_impl< T, tag::addressing_index<2> >::result_type type ;
  } ;
  template< typename T >
  struct size_minor {
      typedef typename detail::size_impl< T, typename addressing_index_minor<T>::type >::result_type type ;
  } ;
  template< typename T >
  struct size_major {
      typedef typename detail::size_impl< T, typename addressing_index_major<T>::type >::result_type type ;
  } ;
}

template< typename T >
typename result_of::size1<T>::type size1( T& t) {
  return detail::size_impl< T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::size1<const T>::type size1( const T& t) {
  return detail::size_impl< const T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::size2<T>::type size2( T& t) {
  return detail::size_impl< T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::size2<const T>::type size2( const T& t) {
  return detail::size_impl< const T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::size3<T>::type size3( T& t) {
  return detail::size_impl< T, tag::addressing_index<3> >::invoke( t ) ;
}
template< typename T >
typename result_of::size3<const T>::type size3( const T& t) {
  return detail::size_impl< const T, tag::addressing_index<3> >::invoke( t ) ;
}
template< typename T >
typename result_of::size_row<T>::type size_row( T& t) {
  return detail::size_impl< T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::size_row<const T>::type size_row( const T& t) {
  return detail::size_impl< const T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::size_column<T>::type size_column( T& t) {
  return detail::size_impl< T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::size_column<const T>::type size_column( const T& t) {
  return detail::size_impl< const T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::size_minor<T>::type size_minor( T& t) {
  return detail::size_impl< T, typename addressing_index_minor<T>::type >::invoke( t ) ;
}
template< typename T >
typename result_of::size_minor<const T>::type size_minor( const T& t) {
  return detail::size_impl< const T, typename addressing_index_minor<T>::type >::invoke( t ) ;
}
template< typename T >
typename result_of::size_minor<T>::type size_major( T& t) {
  return detail::size_impl< T, typename addressing_index_major<T>::type >::invoke( t ) ;
}
template< typename T >
typename result_of::size_minor<const T>::type size_major( const T& t) {
  return detail::size_impl< const T, typename addressing_index_major<T>::type >::invoke( t ) ;
}

/*
#define GENERATE_SIZE_INDEX( z, which, unused ) \
GENERATE_FUNCTIONS( size, which, tag::addressing_index<which> )

BOOST_PP_REPEAT_FROM_TO(1,3,GENERATE_SIZE_INDEX,~)

GENERATE_FUNCTIONS( size, _row, tag::addressing_index<1> )
GENERATE_FUNCTIONS( size, _column, tag::addressing_index<2> )
GENERATE_FUNCTIONS( size, _major, typename addressing_index_major<T>::type )
GENERATE_FUNCTIONS( size, _minor, typename addressing_index_minor<T>::type )
  */

//
// Overloads for free template functions size_row( x, tag ), 
// Here, tag is assumed to be either one of
// tag::transpose, tag::no_transpose, or tag::conjugate
//
namespace result_of {

template< typename T, typename TransTag >
struct size_row_op {
    typedef typename size<
        T,
        typename addressing_index_trans< tag::addressing_index<1>, TransTag >::type
    >::type type;
};

template< typename T, typename TransTag >
struct size_column_op {
    typedef typename size< T, 
        typename addressing_index_trans< tag::addressing_index<2>, TransTag >::type >::type type;
};

} // namespace result_of

template< typename T, typename Tag >
inline typename result_of::size_row_op< const T, Tag >::type
size_row_op( const T& t, Tag ) {
    return bindings::size( t, typename addressing_index_trans< tag::addressing_index<1>, Tag >::type() );
}

template< typename T, typename Tag >
inline typename result_of::size_row_op< const T, Tag >::type
size_column_op( const T& t, Tag ) {
    return bindings::size( t, typename addressing_index_trans< tag::addressing_index<2>, Tag >::type() );
}

} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
