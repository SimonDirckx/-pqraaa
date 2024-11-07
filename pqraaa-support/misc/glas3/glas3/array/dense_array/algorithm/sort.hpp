//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_sort_hpp
#define glas3_array_dense_array_algorithm_sort_hpp

#include <glas3/array/dense_array/type/iterator.hpp>
#include <iterator>
#include <type_traits>

//namespace glas3 { namespace sort_detail {
//
//  template <typename E, typename I>
//  class sort_iterator {
//    public:
//      sort_iterator( E const& e, I const& ind, std::ptrdiff_t pos )
//      : e_( boost::make_shared< E >( std::move( e.shallow_copy() ) ) )
//      , ind_( boost::make_shared< I >( std::move( ind.shallow_copy() ) ) )
//      , pos_( pos )
//      {}
//
//    public:
//      struct value_type {
//        typedef typename E::value_type first_type ;
//        typedef typename I::value_type second_type ;
//
//        value_type( first_type e, second_type ind )
//        : e_( e )
//        , ind_( ind )
//        {}
//
//        bool operator<( value_type const& that ) const { return e_ < that.e_ ; }
//        bool operator>( value_type const& that ) const { return e_ > that.e_ ; }
//        bool operator==( value_type const& that ) const { return e_ == that.e_ ; }
//        bool operator!=( value_type const& that ) const { return e_ != that.e_ ; }
//
//        first_type  e_ ;
//        second_type ind_ ;
//      } ;
//
//      typedef value_type const_reference ;
//
//      struct reference {
//        typedef typename E::value_type& first_type ;
//        typedef typename I::value_type& second_type ;
//
//        reference( first_type e, second_type ind )
//        : e_( e )
//        , ind_( ind )
//        {}
//
////        public:
////	    // Copy constructor
////        reference ( reference const& that ) = default ;
//
//        // Move constructor
//        reference ( reference && that ) = default;
//
//        public:
//        // Move assignment
//        reference& operator= ( reference&& that ) = default;
//
//        reference& operator=( value_type const& that ) { e_ = that.e_ ; ind_ = that.ind_ ; return *this ; }
//        reference& operator=( reference const& that ) { e_ = that.e_ ; ind_ = that.ind_ ; return *this ; }
//
//        operator value_type() const { return value_type( e_, ind_ ) ; }
//
//        bool operator<( reference const& that ) const { return e_ < that.e_ ; }
//        bool operator>( reference const& that ) const { return e_ > that.e_ ; }
//        bool operator==( reference const& that ) const { return e_ == that.e_ ; }
//        bool operator!=( reference const& that ) const { return e_ != that.e_ ; }
//
//        bool operator<( const_reference const& that ) const { return e_ < that.e_ ; }
//        bool operator>( const_reference const& that ) const { return e_ > that.e_ ; }
//        bool operator==( const_reference const& that ) const { return e_ == that.e_ ; }
//        bool operator!=( const_reference const& that ) const { return e_ != that.e_ ; }
//
//        first_type  e_ ;
//        second_type ind_ ;
//      } ;
//      typedef value_type*                                        pointer ;
//
//      typedef std::random_access_iterator_tag                    iterator_category ;
//      typedef std::ptrdiff_t                                     difference_type ;
//
//    public:
//      sort_iterator( sort_iterator const& that )
//      : e_( that.e_ )
//      , ind_( that.ind_ )
//      , pos_( that.pos_ )
//      {}
//
//      sort_iterator& operator=( sort_iterator const& that ) { pos_ = that.pos_ ; return *this ; }
//
//    public:
//      const_reference operator*() const { return const_reference( (*e_)[pos_], (*ind_)[pos_] ) ; }
//      //reference operator*() { return reference( glas3::at( e_,  pos_ ), glas3::at( ind_, pos_ ) ) ; }
//      reference operator*() { return reference( (*e_)[pos_], (*ind_)[pos_] ) ; }
//
//      const_reference operator[]( std::ptrdiff_t i ) const { return const_reference( (*e_)[pos_ + i], (*ind_)[pos_ + i] ) ; }
//      reference operator[]( std::ptrdiff_t i ) { return reference( (*e_)[pos_ + i], (*ind_)[pos_ + i] ) ; }
//
//      sort_iterator& operator++() { ++pos_ ; return *this ; }
//      sort_iterator& operator--() { --pos_ ; return *this ; }
//
//      sort_iterator operator++( int ) { sort_iterator that( *this ) ; ++pos_ ; return that ; }
//      sort_iterator operator--( int ) { sort_iterator that( *this ) ; --pos_ ; return that ; }
//
//      sort_iterator& operator+=( std::ptrdiff_t i ) { pos_ += i ; return *this ; }
//      sort_iterator& operator-=( std::ptrdiff_t i ) { pos_ -= i ; return *this ; }
//      difference_type operator-( sort_iterator const& that ) const { return pos_ - that.pos_ ; }
//      sort_iterator operator-( difference_type i ) const { return sort_iterator( *e_, *ind_, pos_ - i ) ; }
//      sort_iterator operator+( difference_type i ) const { return sort_iterator( *e_, *ind_, pos_ + i ) ; }
//
//      bool operator<( sort_iterator const& that ) const { return pos_ < that.pos_ ; }
//      bool operator>( sort_iterator const& that ) const { return pos_ > that.pos_ ; }
//      bool operator==( sort_iterator const& that ) const { return pos_ == that.pos_ ; }
//      bool operator!=( sort_iterator const& that ) const { return pos_ != that.pos_ ; }
//
//      friend sort_iterator operator+( difference_type i, sort_iterator const& s ) { return sort_iterator( *(s.e_), *(s.ind_), i + s.pos_ ) ; }
//      friend sort_iterator operator-( difference_type i, sort_iterator const& s ) { return sort_iterator( *(s.e_), *(s.ind_), i - s.pos_ ) ; }
//
//    private:
//      boost::shared_ptr< E >   e_ ;
//      boost::shared_ptr< I >   ind_ ;
//      std::ptrdiff_t           pos_ ;
//  } ;
//
//} } // namespace glas3::sort_detail


namespace glas3 {

//
// Sort (default: ascending order): undefined behavior for indirect DenseArray's with recurring entries; DenseArray should be assignable
//

template <typename E>
typename std::enable_if< is< DenseArray, E >::value >::type sort( E const& e ) {
	std::sort( glas3::iterator<E>(e, 0), glas3::iterator<E>( e, e.size() ) ) ;
}

template <typename E, typename Compare>
typename std::enable_if< is< DenseArray, E >::value >::type sort( E const& e, Compare const& c) {
	std::sort( glas3::iterator<E>(e, 0), glas3::iterator<E>( e, e.size() ), c ) ;
}

template <typename E>
typename std::enable_if< is< DenseArray, E >::value >::type sort_ascending( E const& e) {
	std::sort( glas3::iterator<E>(e, 0), glas3::iterator<E>( e, e.size() ) ) ;
}

template <typename E>
typename std::enable_if< is< DenseArray, E >::value >::type sort_descending( E const& e) {
	std::sort( glas3::iterator<E>(e, 0), glas3::iterator<E>( e, e.size() ), std::function< bool ( typename E::value_type const&, typename E::value_type const& ) >( [] ( typename E::value_type const& x, typename E::value_type const& y ) { return x > y ; } ) ) ;
}

//
// Sort in ascending order and let ind follow: COMPILE ERROR??
//

//  template <typename E, typename I>
//  void sort( E const& e, I const& ind ) {
//    assert( e.size() == ind.size() ) ;
//    std::sort( glas3::sort_detail::sort_iterator<E, I>( e, ind, 0 ), glas3::sort_detail::sort_iterator<E, I>( e, ind, e.size() ) ) ;
//  }

} // namespace glas2

#endif
