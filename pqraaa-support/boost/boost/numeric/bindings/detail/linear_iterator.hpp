//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_DETAIL_LINEAR_ITERATOR_HPP
#define BOOST_NUMERIC_BINDINGS_DETAIL_LINEAR_ITERATOR_HPP

#include <iterator>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename StrideType >
class linear_iterator {
  public:
    typedef std::random_access_iterator_tag category ;
    typedef std::integral_constant<int, StrideType::value> stride_type ;
    typedef T* pointer ;

    linear_iterator()
    : it_( 0 )
    {}

    linear_iterator( pointer p, StrideType ignore )
    : it_( p )
    {}

  public:
    linear_iterator operator++() { it_ += stride_type::value ; return *this ; }
    linear_iterator operator+=( int n ) { it_ += n*stride_type::value ; return *this ; }

    linear_iterator operator--() { it_ -= stride_type::value ; return *this ; }
    linear_iterator operator-=( int n ) { it_ -= n*stride_type::value ; return *this ; }

    pointer operator*() const { return it_ ; }

  private:
    pointer it_ ;
};

template< typename T >
class linear_iterator< T, std::ptrdiff_t > {
  public:
    typedef std::random_access_iterator_tag category ;
    typedef T* pointer ;

    linear_iterator()
    : it_( 0 )
    , m_stride( 0 )
    {}

    linear_iterator( pointer p, std::ptrdiff_t stride )
    : it_( p )
    , m_stride( stride )
    {}

  public:
    linear_iterator operator++() { it_ += m_stride ; return *this ; }
    linear_iterator operator+=( int n ) { it_ += n*m_stride ; return *this ; }

    linear_iterator operator--() { it_ -= m_stride ; return *this ; }
    linear_iterator operator-=( int n ) { it_ -= n*m_stride ; return *this ; }

    pointer operator*() const { return it_ ; }

  private:
    pointer it_ ;
    std::ptrdiff_t m_stride;
};

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
