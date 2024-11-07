//  (C) Copyright Karl Meerbergen 2022.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_diagonal_matrix_hpp
#define glas2_matrix_type_diagonal_matrix_hpp

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/algorithm/select.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename V>
  class diagonal_matrix {
    public:
      typedef typename V::value_type value_type ;
      typedef typename V::size_type  size_type ;

    public:
      //diagonal_matrix( V const& v, size_type diag=0, size_type n_rows=v.size(), size_type n_cols=v.size() )
      diagonal_matrix( V const& v )
      : v_(v)
      , diag_( 0 )
      , n_rows_( v.size() )
      , n_cols_( v.size() )
      {
        assert(diag_==0) ;
        assert(n_rows_=v_.size()) ;
        assert(n_cols_=v_.size()) ;
        assert( diag_+v_.size() >= n_cols_ ) ;
        assert( -diag_+v_.size() >= n_rows_ ) ;
      }

      // Copy reference !!
      diagonal_matrix( diagonal_matrix const& that )
      : v_( that.v_ )
      , diag_( that.diag_ )
      , n_rows_( that.n_rows_ )
      , n_cols_( that.n_cols_ )
      {}

    public:
      // Matrix
      size_type num_rows() const { return n_rows_ ; }
      size_type num_columns() const { return n_cols_ ; }

      V const& diagonal() const { return v_ ; }

      value_type operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows() ) ;
        assert( j>=0 && j<num_columns() ) ;
        if (i-diag_==j) return v_(i) ;
        else return 0 ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< diagonal_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) {
        return matrix_selection< diagonal_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      V         v_ ;
      size_type diag_ ;
      size_type n_rows_ ;
      size_type n_cols_ ;
  } ;


  template <typename V>
  struct glas_concept< diagonal_matrix<V> >
  : Matrix
  {};

} // namespace glas2


#endif
