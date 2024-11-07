//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_algorithm_sort_hpp
#define glas2_vector_algorithm_sort_hpp

#include <glas2/vector/type/iterator.hpp>
#include <iterator>
#include <type_traits>

namespace glas2 { namespace sort_detail {

  template <typename E, typename I>
  class sort_iterator {
    public:
      sort_iterator( E& e, I& ind, std::ptrdiff_t pos )
      : e_( e )
      , ind_( ind )
      , pos_( pos )
      {}

    public:
      struct value_type {
        typedef typename E::value_type first_type ;
        typedef typename I::value_type second_type ;

        value_type( first_type e, second_type ind )
        : e_( e )
        , ind_( ind )
        {}

        bool operator<( value_type const& that ) const { return e_ < that.e_ ; }
        bool operator>( value_type const& that ) const { return e_ > that.e_ ; }
        bool operator>=( value_type const& that ) const { return e_ >= that.e_ ; }
        bool operator==( value_type const& that ) const { return e_ == that.e_ ; }
        bool operator!=( value_type const& that ) const { return e_ != that.e_ ; }

        first_type  e_ ;
        second_type ind_ ;
      } ;

      typedef value_type const_reference ;

      struct reference {
        typedef typename E::value_type& first_type ;
        typedef typename I::value_type& second_type ;

        reference( first_type e, second_type ind )
        : e_( e )
        , ind_( ind )
        {}

        reference& operator=( value_type const& that ) { e_ = that.e_ ; ind_ = that.ind_ ; return *this ; }
        reference& operator=( reference const& that ) { e_ = that.e_ ; ind_ = that.ind_ ; return *this ; }

        operator value_type() const { return value_type( e_, ind_ ) ; }

        bool operator<( reference const& that ) const { return e_ < that.e_ ; }
        bool operator>( reference const& that ) const { return e_ > that.e_ ; }
        bool operator>=( reference const& that ) const { return e_ >= that.e_ ; }
        bool operator==( reference const& that ) const { return e_ == that.e_ ; }
        bool operator!=( reference const& that ) const { return e_ != that.e_ ; }

        bool operator<( const_reference const& that ) const { return e_ < that.e_ ; }
        bool operator>( const_reference const& that ) const { return e_ > that.e_ ; }
        bool operator>=( const_reference const& that ) const { return e_ >= that.e_ ; }
        bool operator==( const_reference const& that ) const { return e_ == that.e_ ; }
        bool operator!=( const_reference const& that ) const { return e_ != that.e_ ; }

        friend void swap( reference x, reference y ) {
          value_type temp = x ;
          x = y; y = temp ;
        }

        first_type  e_ ;
        second_type ind_ ;
      } ;
      typedef value_type*                                        pointer ;

      typedef std::random_access_iterator_tag                    iterator_category ;
      typedef std::ptrdiff_t                                     difference_type ;

    public:
      sort_iterator( sort_iterator const& that ) 
      : e_( that.e_ )
      , ind_( that.ind_ )
      , pos_( that.pos_ )
      {}


      sort_iterator& operator=( sort_iterator const& that ) { pos_ = that.pos_ ; return *this ; }

    public:
      const_reference operator*() const { return cont_reference( e_(pos_), ind_(pos_) ) ; }
      //reference operator*() { return reference( glas::at( e_,  pos_ ), glas::at( ind_, pos_ ) ) ; }
      reference operator*() {
        assert( pos_>=0 && pos_<e_.size() ) ;
        return reference( e_(pos_), ind_(pos_) ) ;
      }

      const_reference operator[]( std::ptrdiff_t i ) const { return const_reference( e_(pos_+i), ind_(pos_+i) ) ; }
      reference operator[]( std::ptrdiff_t i ) { return reference( e_(pos_+i), ind_(pos_+i) ) ; }

      sort_iterator& operator++() { ++pos_ ; return *this ; }
      sort_iterator& operator--() { --pos_ ; return *this ; }

      sort_iterator operator++( int ) { sort_iterator that( *this ) ; ++pos_ ; return that ; }
      sort_iterator operator--( int ) { sort_iterator that( *this ) ; --pos_ ; return that ; }

      sort_iterator& operator+=( std::ptrdiff_t i ) { pos_ += i ; return *this ; }
      sort_iterator& operator-=( std::ptrdiff_t i ) { pos_ -= i ; return *this ; }
      difference_type operator-( sort_iterator const& that ) const { return pos_ - that.pos_ ; }
      sort_iterator operator-( difference_type i ) const { return sort_iterator( e_, ind_, pos_-i ) ; }
      sort_iterator operator+( difference_type i ) const { return sort_iterator( e_, ind_, pos_+i ) ; }

      bool operator<( sort_iterator const& that ) const { return pos_ < that.pos_ ; }
      bool operator>( sort_iterator const& that ) const { return pos_ > that.pos_ ; }
      bool operator>=( sort_iterator const& that ) const { return pos_ >= that.pos_ ; }
      bool operator==( sort_iterator const& that ) const { return pos_ == that.pos_ ; }
      bool operator!=( sort_iterator const& that ) const { return pos_ != that.pos_ ; }

      friend sort_iterator operator+( difference_type i, sort_iterator const& s ) { return sort_iterator( s.e_, s.ind_, i+s.pos_ ) ; }
      friend sort_iterator operator-( difference_type i, sort_iterator const& s ) { return sort_iterator( s.e_, s.ind_, i-s.pos_ ) ; }

    private:
      E&             e_ ;
      I&             ind_ ;
      std::ptrdiff_t pos_ ;
  } ;

} } // namespace glas2::sort_detail


namespace glas2 {

  //
  // Sort in ascending order
  //

  template <typename E>
  E sort( E e ) {
    std::sort( glas2::iterator<E>(e,0), glas2::iterator<E>( e, e.size() ) ) ;
    return e ;
  }

  //
  // Sort with binary comparison function
  //
 
  template <typename E, typename C>
  E sort_comp( E e, C comp ) {
    std::sort( glas2::iterator<E>(e,0), glas2::iterator<E>( e, e.size() ), comp) ;
    return e ;
  }

  //
  // Sort in ascending order and let ind follow
  //

  template <typename E, typename I>
  E sort( E e, I ind ) {
    assert( e.size()==ind.size() ) ;
    std::sort( glas2::sort_detail::sort_iterator<E,I>( e, ind, 0 ), glas2::sort_detail::sort_iterator<E,I>( e, ind, e.size() ) ) ;
    return e ;
  }

} // namespace glas2

#endif
