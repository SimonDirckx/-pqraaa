//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_STRIDE_HPP
#define BOOST_NUMERIC_BINDINGS_STRIDE_HPP

#include <boost/numeric/bindings/size.hpp>
#include <type_traits>
#include <typeinfo>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename AddressingIndex, typename Enable = void >
struct stride_impl {

    typedef tag::stride_type< AddressingIndex::value > key_type;
    typedef typename result_of_get< T, key_type >::type result_type;

    static result_type invoke( const T& t ) {
        return get< key_type >( t );
    }

};

//
// Strides for ranks outside the scope of the object are fixed at
// the dot product of its existing sizes and strides.
// 
// Object    rank    result of strideN(), with N > rank
// scalar    0       1
// vector    1       size1 * stride1
// matrix    2       size1 * stride1 + size2 * stride2
// tensor    N       sum_i( size_i, stride_i )  (dot( size, stride))
//
// Iff size_i and stride_i are integral constants, results will be known at
// compile time. Otherwise, the result_type will be std::ptrdiff_t.
//
template< typename T, typename AddressingIndex >
struct stride_impl< T, AddressingIndex,
        typename std::enable_if<
            std::is_same< rank<T>, tag::scalar >::value
        >::type > {

    typedef typename std::integral_constant<int,1> result_type;

    static result_type invoke( const T& t ) {
        return result_type();
    }

};


template< typename T, typename State, typename AddressingIndex >
struct fold_stride_size {

    typedef tag::addressing_index< AddressingIndex::value > index_type;
    typedef typename result_of::size< T, index_type >::type size_type;
    typedef typename stride_impl< T, index_type >::result_type stride_type;

    typedef typename std::conditional<
                (std::is_same< State, std::ptrdiff_t >::value || std::is_same< size_type, std::ptrdiff_t >::value)
                || std::is_same< stride_type, std::ptrdiff_t >::value,
                std::ptrdiff_t,
                decltype( State() + size_type()*stride_type() )
                >::type type ;
//    typedef typename mpl::if_<
//        mpl::or_<
//            is_same< State, std::ptrdiff_t >,
//            is_same< size_type, std::ptrdiff_t >,
//            is_same< stride_type, std::ptrdiff_t >
//        >,
//        std::ptrdiff_t,
//        mpl::plus<
//            State, 
//            mpl::times< 
//                size_type,
//                stride_type
//            >
//        >
//    >::type type;

};

//
// If Result isn't a ptrdiff_t, just invoke the integral constant
// and return that. Otherwise, runtime stuff is involved, so we'll
// have to evaluate sum_i( size_i, stride_i ).
//
template< typename T, typename Result, int AddressingIndex >
struct apply_fold {
    static Result invoke( const T& t ) {
        return Result();
    }
};

template< typename T, int AddressingIndex >
struct apply_fold< T, std::ptrdiff_t, AddressingIndex > {

    static std::ptrdiff_t invoke( const T& t ) {
        return size( t, tag::addressing_index< AddressingIndex >() ) * 
            stride_impl< T, tag::addressing_index< AddressingIndex > >::invoke( t ) +
            apply_fold< T, std::ptrdiff_t, AddressingIndex-1 >::invoke( t );
    }

};

template< typename T >
struct apply_fold< T, std::ptrdiff_t, 0 > {

    static std::ptrdiff_t invoke( const T& ) {
        return 0;
    }

};


template <typename T, int AddressingIndex >
struct stride_sum_impl {
  typedef tag::addressing_index< AddressingIndex > index_type;
  typedef typename result_of::size< T, index_type >::type size_type;
  typedef typename stride_impl< T, index_type >::result_type stride_type;

  typedef decltype( typename stride_sum_impl<T, AddressingIndex-1>::result_type() + size_type()*stride_type() ) result_type ;

  static result_type invoke( const T& t ) {
    return stride_sum_impl<T, AddressingIndex-1>::invoke( t ) + stride_impl< T, tag::addressing_index<AddressingIndex> >::invoke( t ) * size( t, tag::addressing_index<AddressingIndex>() ) ;
  }
} ;

template <typename T >
struct stride_sum_impl<T, 1>
{
  typedef tag::addressing_index< 1 > index_type ;
  typedef typename result_of::size< T, index_type >::type size_type;
  typedef typename stride_impl< T, index_type >::result_type stride_type;

  typedef decltype( size_type()*stride_type() ) result_type ;

  static result_type invoke( const T& t ) {
    return stride_impl< T, index_type>::invoke( t ) * size( t, index_type() ) ;
  }
} ;

template <typename T >
struct stride_sum_impl<T, 0>
{
  typedef tag::addressing_index< 0 > index_type ;
  typedef typename result_of::size< T, index_type >::type size_type;
  typedef typename stride_impl< T, index_type >::result_type stride_type;

  typedef decltype( size_type()*stride_type() ) result_type ;

  static result_type invoke( const T& t ) {
    return stride_impl< T, index_type>::invoke( t ) * size( t, index_type() ) ;
  }
} ;

// Could be made generic for dimensions > 2,
//  but not enough time right now


template< typename T, typename AddressingIndex >
struct stride_impl< T, AddressingIndex,
        typename std::enable_if< (rank<T>::value > tag::scalar::value) && (AddressingIndex::value > rank<T>::value) >::type > {
//        typename boost::enable_if<
//            mpl::and_<
//                 mpl::greater< rank<T>, tag::scalar >,
//                 mpl::greater< AddressingIndex, rank<T> >
//           >
//        >::type > {

    typedef typename stride_sum_impl< T, rank<T>::value >::result_type result_type ;
/*
    typedef mpl::range_c< int, 1, rank<T>::value+1 > index_range;
    typedef typename mpl::fold<
        index_range,
        mpl::int_< 0 >,
        fold_stride_size<
            T,
            mpl::_1,
            mpl::_2
        >
    >::type result_type;

    static result_type invoke( const T& t ) {
        return apply_fold< T, result_type, rank<T>::value >::invoke( t );
    }
*/
    static result_type invoke( const T& t ) {
        return stride_sum_impl< T, rank<T>::value >::invoke( t );
    }

};

} // namespace detail

namespace result_of {

template< typename T, typename Tag = tag::addressing_index<1> >
struct stride {
    static_assert( (is_tag<Tag>::value) );
    typedef typename detail::stride_impl< T, Tag >::result_type type;
};

} // namespace result_of


//
// Overloads for free template functions stride( x, tag ), 
//
template< typename T, typename Tag >
inline typename result_of::stride< const T, Tag >::type
stride( const T& t, Tag ) {
    return detail::stride_impl< const T, Tag >::invoke( t );
}

// Overloads for free template functions stride( x )
// Valid for types with rank <= 1 (scalars, vectors)
// In theory, we could provide overloads for matrices here, too, 
// if their minimal_rank is at most 1.

template< typename T >
inline typename
//boost::enable_if<
//    mpl::less< rank<T>, mpl::int_<2> >,
std::enable_if<
    rank<T>::value < 2,
    typename result_of::stride< const T, tag::addressing_index<1> >::type
>::type
stride( const T& t ) {
    return detail::stride_impl<const T, tag::addressing_index<1> >::invoke( t );
}

namespace result_of {
  template< typename T >
  struct stride1 {
      typedef typename detail::size_impl< T, tag::addressing_index<1> >::result_type type ;
  } ;
  template< typename T >
  struct stride2 {
      typedef typename detail::size_impl< T, tag::addressing_index<2> >::result_type type ;
  } ;
  template< typename T >
  struct stride3 {
      typedef typename detail::size_impl< T, tag::addressing_index<3> >::result_type type ;
  } ;
  template< typename T >
  struct stride_row {
      typedef typename detail::size_impl< T, tag::addressing_index<1> >::result_type type ;
  } ;
  template< typename T >
  struct stride_column {
      typedef typename detail::size_impl< T, tag::addressing_index<2> >::result_type type ;
  } ;
  template< typename T >
  struct stride_major {
      typedef typename detail::size_impl< T, typename addressing_index_major<T>::type >::result_type type ;
  } ;
  template< typename T >
  struct stride_minor {
      typedef typename detail::size_impl< T, typename addressing_index_minor<T>::type >::result_type type ;
  } ;
}

template< typename T >
typename result_of::stride1<T>::type stride1( T& t) {
  return detail::stride_impl< T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride1<const T>::type stride1( const T& t) {
  return detail::stride_impl< const T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride2<T>::type stride2( T& t) {
  return detail::stride_impl< T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride2<const T>::type stride2( const T& t) {
  return detail::stride_impl< const T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride3<T>::type stride3( T& t) {
  return detail::stride_impl< T, tag::addressing_index<3> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_row<T>::type stride_row( T& t) {
  return detail::stride_impl< T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_row<const T>::type stride_row( const T& t) {
  return detail::stride_impl< const T, tag::addressing_index<1> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_column<T>::type stride_column( T& t) {
  return detail::stride_impl< T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_column<const T>::type stride_column( const T& t) {
  return detail::stride_impl< const T, tag::addressing_index<2> >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_major<T>::type stride_major( T& t) {
  return detail::stride_impl< T, typename addressing_index_major<T>::type >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_major<const T>::type stride_major( const T& t) {
  return detail::stride_impl< const T, typename addressing_index_major<T>::type >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_minor<T>::type stride_minor( T& t) {
  return detail::stride_impl< T, typename addressing_index_minor<T>::type >::invoke( t ) ;
}
template< typename T >
typename result_of::stride_minor<const T>::type stride_minor( const T& t) {
  return detail::stride_impl< const T, typename addressing_index_minor<T>::type >::invoke( t ) ;
}

/*
#define GENERATE_STRIDE_INDEX( z, which, unused ) \
GENERATE_FUNCTIONS( stride, which, tag::addressing_index<which> )

BOOST_PP_REPEAT_FROM_TO(1,3,GENERATE_STRIDE_INDEX,~)

GENERATE_FUNCTIONS( stride, _row, tag::addressing_index<1> )
GENERATE_FUNCTIONS( stride, _column, tag::addressing_index<2> )
GENERATE_FUNCTIONS( stride, _major, typename addressing_index_major<T>::type )
GENERATE_FUNCTIONS( stride, _minor, typename addressing_index_minor<T>::type )
*/

} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
