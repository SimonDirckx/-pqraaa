//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_iterator_hpp
#define glas2_vector_type_iterator_hpp

namespace glas2 {

  template <typename V>
  class iterator {
    public:
      typedef typename V::value_type value_type ;
      typedef typename V::size_type  size_type ;
      typedef value_type const&      const_reference ;
      typedef value_type&            reference ;
      typedef value_type*            pointer ;

      typedef std::random_access_iterator_tag  iterator_category ;
      typedef std::ptrdiff_t                   difference_type ;

    public:
      inline iterator( V v, size_type i )
      : v_( v )
      , index_( i )
      {}

    public:
      iterator& operator=( iterator const& that ) {
        index_ = that.index_ ;
        return *this ;
      }

    public:
      const_reference operator*() const { return v_(index_) ; }
      reference operator*() { return v_(index_) ; }

      const_reference operator[]( std::ptrdiff_t i ) const { return v_(index_+i) ; }
      reference operator[]( std::ptrdiff_t i ) { return v_(index_+i) ; }

      iterator& operator++() { ++index_ ; return *this ; }
      iterator& operator--() { --index_ ; return *this ; }

      iterator operator++( int ) { iterator that( *this ) ; ++index_ ; return that ; }
      iterator operator--( int ) { iterator that( *this ) ; --index_ ; return that ; }

      iterator& operator+=( std::ptrdiff_t i ) { index_ += i ; return *this ; }
      iterator& operator-=( std::ptrdiff_t i ) { index_ -= i ; return *this ; }
      difference_type operator-( iterator const& that ) const { return index_ - that.index_ ; }
      iterator operator-( difference_type i ) const { return iterator( v_, index_-i ) ; }

      template<typename I>
      iterator operator+( I i ) const { return iterator( v_, index_+i ) ; }

      bool operator<( iterator const& that ) const { return index_ < that.index_ ; }
      bool operator<=( iterator const& that ) const { return index_ <= that.index_ ; }
      bool operator>=( iterator const& that ) const { return index_ >= that.index_ ; }
      bool operator>( iterator const& that ) const { return index_ > that.index_ ; }
      bool operator==( iterator const& that ) const { return index_ == that.index_ ; }
      bool operator!=( iterator const& that ) const { return index_ != that.index_ ; }

      friend iterator operator+( difference_type i, iterator const& s ) { return iterator( s.v_, i+s.index_ ) ; }
      friend iterator operator-( difference_type i, iterator const& s ) { return iterator( s.v_, i-s.index_ ) ; }

    private:
      V         v_ ;
      size_type index_ ;
  } ;

} // namespace glas2

#endif
