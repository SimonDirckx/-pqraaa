//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_algorithm_sort_and_compress_hpp
#define glas2_sparse_algorithm_sort_and_compress_hpp

#include <iterator>
#include <algorithm>

namespace glas2 {

  namespace detail {
    template <typename T>
    struct less {
      bool operator() ( T const& t1, T const& t2 ) const {
        return t1.row_ < t2.row_ || (t1.row_==t2.row_ && t1.column_ < t2.column_ ) ;
      }
    } ;

    template <typename COO>
    //template <typename COO, template <typename T> typename SmallerThan=less>
    class coo_iterator {
      private:
        typedef typename COO::index_type index_type ;

      public:
        typedef typename COO::size_type size_type ;
        typedef size_type difference_type ;

        typedef std::random_access_iterator_tag iterator_category ;

        struct value_type {
          typename COO::value_type value_ ;
          typename COO::index_type row_ ;
          typename COO::index_type column_ ;

          value_type( typename COO::value_type const& value, typename COO::index_type const& row, typename COO::index_type const& column )
          : value_( value )
          , row_( row )
          , column_( column )
          {}

          bool operator<( value_type const& that ) const {
            return row_ < that.row_ || (row_==that.row_ && column_ < that.column_ ) ;
          }

          bool operator<=( value_type const& that ) const {
            return row_ < that.row_ || (row_==that.row_ && column_ <= that.column_ ) ;
          }

          bool operator==( value_type const& that ) const {
            return row_==that.row_ && column_ == that.column_ ;
          }

          bool operator!=( value_type const& that ) const {
            return row_!=that.row_ || column_ != that.column_ ;
          }

        } ;

        typedef value_type  const_reference ;
        typedef value_type* pointer ;

        struct reference {
          reference( typename COO::value_type& value, typename COO::index_type& row, typename COO::index_type& column )
          : value_( value )
          , row_( row )
          , column_( column )
          {}

          reference( reference const& that )
          : value_( that.value_ )
          , row_( that.row_ )
          , column_( that.column_ )
          {}

          reference& operator=( reference const& that ) {
            value_ = that.value_ ;
            row_ = that.row_ ;
            column_ = that.column_ ;
            return *this ;
          }

          reference& operator=( value_type const& that ) {
            value_ = that.value_ ;
            row_ = that.row_ ;
            column_ = that.column_ ;
            return *this ;
          }

          operator value_type() const { return value_type( value_, row_, column_ ) ; }

          friend void swap( reference one, reference two ) {
            typename COO::value_type v = one.value_ ; one.value_ = two.value_ ; two.value_ = v ;
            typename COO::index_type i = one.row_ ; one.row_ = two.row_ ; two.row_ = i ;
            i = one.column_ ; one.column_ = two.column_ ; two.column_ = i ;
          }

          bool operator<( reference const& that ) const {
            return row_ < that.row_ || (row_==that.row_ && column_ < that.column_ ) ;
          }
          bool operator<( value_type const& that ) const {
            return row_ < that.row_ || (row_==that.row_ && column_ < that.column_ ) ;
          }

          bool operator<=( reference const& that ) const {
            return row_ < that.row_ || (row_==that.row_ && column_ <= that.column_ ) ;
          }
          bool operator<=( value_type const& that ) const {
            return row_ < that.row_ || (row_==that.row_ && column_ <= that.column_ ) ;
          }

          bool operator==( reference const& that ) const {
            return row_==that.row_ && column_ == that.column_ ;
          }
          bool operator==( value_type const& that ) const {
            return row_==that.row_ && column_ == that.column_ ;
          }

          bool operator!=( reference const& that ) const {
            return row_!=that.row_ || column_ != that.column_ ;
          }
          bool operator!=( value_type const& that ) const {
            return row_!=that.row_ || column_ != that.column_ ;
          }

          typename COO::value_type& value_ ;
          typename COO::index_type& row_ ;
          typename COO::index_type& column_ ;
        } ;

      public:
        coo_iterator( COO& coo, size_type i )
        : coo_( coo )
        , i_( i )
        {}

        reference operator*() {
          return reference( coo_.data()(i_), coo_.row_indices()(i_), coo_.column_indices()(i_) ) ;
        }

        const_reference operator*() const{
          value_type v;
          v.value_ = coo_.data()(i_) ;
          v.row_ = coo_.row(i_) ;
          v.column_ = coo_.column(i_) ;
          return v ;
        }

      public:
        coo_iterator( coo_iterator const& that )
        : coo_( that.coo_ )
        , i_( that.i_ )
        {}

        coo_iterator& operator=( coo_iterator const& that ) {
          i_ = that.i_ ;
          return *this ;
        }

      public:
        coo_iterator swap( coo_iterator& that ) {
          typename COO::value_type d = coo_.data()( i_ ) ;
          coo_.data()( i_ ) = that.coo_.data()( that.i_ ) ;
          that.coo_.data()( that.i_ ) = d ;

          index_type r = coo_.row( i_ ) ;
          coo_.row()( i_ ) = that.coo_.row()( that.i_ ) ;
          that.coo_.row()( that.i_ ) = r ;

          index_type c = coo_.column( i_ ) ;
          coo_.column()( i_ ) = that.coo_.column()( that.i_ ) ;
          that.coo_.column()( that.i_ ) = c ;
        }

      public:
        bool operator==( coo_iterator const& that ) const {
          return i_ == that.i_ ;
        }

        bool operator!=( coo_iterator const& that ) const {
          return i_ != that.i_ ;
        }

        bool operator< ( coo_iterator const& that ) const {
          return i_ < that.i_ ;
        }

        bool operator<= ( coo_iterator const& that ) const {
          return i_ <= that.i_ ;
        }

        bool operator>= ( coo_iterator const& that ) const {
          return i_ >= that.i_ ;
        }

        bool operator> ( coo_iterator const& that ) const {
          return i_ > that.i_ ;
        }

      public:
        difference_type operator-( coo_iterator const& that ) const {
          return i_ - that.i_ ;
        }

      public:
        coo_iterator operator+( int n ) {
          coo_iterator that( *this ) ;
          that.i_ += n ;
          return that ;
        }

        coo_iterator operator-( int n ) {
          coo_iterator that( *this ) ;
          that.i_ -= n ;
          return that ;
        }

        coo_iterator& operator+=( int n ) {
          i_ += n ;
          return *this ;
        }

        coo_iterator& operator-=( int n ) {
          i_ -= n ;
          return *this ;
        }

        coo_iterator& operator++() {
          i_ ++ ;
          return *this ;
        }

        coo_iterator& operator--() {
          i_ -- ;
          return *this ;
        }

      private:
        COO&      coo_ ;
        size_type i_ ;
    };
  } // namespace detail

  template <typename COO>
  COO& sort_and_compress( COO& coo ) {
    std::sort( detail::coo_iterator<COO>( coo, 0 ) , detail::coo_iterator<COO>( coo, coo.num_nz() ) ) ;

    // Remove doubles

    return coo ;
  }

}

#endif
