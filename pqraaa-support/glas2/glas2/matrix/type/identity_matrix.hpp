//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_identity_matrix_hpp
#define glas2_matrix_type_identity_matrix_hpp

#include <glas2/matrix/algorithm/assign.hpp>
#include <glas2/matrix/algorithm/select.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/unit_vector.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T=double>
  class identity_matrix {
    public:
      typedef T value_type ;
      typedef int size_type ;
      typedef any_orientation orientation ;

    public:
      identity_matrix()
      : num_rows_(0)
      , num_columns_(0)
      {}

      identity_matrix( size_type m, size_type n )
      : num_rows_( m )
      , num_columns_( n )
      {}

      // Copy reference !!
      identity_matrix( identity_matrix const& that )
      : num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      {}

    public:
      void reset( size_type m, size_type n ) {
        num_rows_ = m ;
        num_columns_ = n ;
      }

    public:
      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }

      value_type operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        if (i==j) return 1 ;
        else return 0 ;
      }

      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< identity_matrix, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< identity_matrix, I1, I2 >::apply( *this, s1, s2 ) ;
      }

    private:
      size_type num_rows_ ;
      size_type num_columns_ ;
  } ;


  template <typename T, typename R, typename I>
  struct matrix_selection< identity_matrix<T>, R, I,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< unit_vector<T>, R>::result_type result_type ;

    static result_type apply( identity_matrix<T> m, R const& r, I col ) {
      assert( col<m.num_columns() && col>=0 ) ;
      return unit_vector<T>( m.num_rows(), col )( r ) ;
    }
  } ;


  template <typename T, typename I, typename R>
  struct matrix_selection< identity_matrix<T>, I, R,
        typename std::enable_if< is<DenseVector,R>::value && std::is_integral<I>::value >::type > {
    typedef typename vector_selection< unit_vector<T>, R>::result_type result_type ;

    static result_type apply( identity_matrix<T> m, I row, R const& r ) {
      assert( row<m.num_rows() && row>=0 ) ;
      return unit_vector<T>( m.num_columns(), row )( r ) ;
    }
  } ;


  template <typename T>
  struct glas_concept< identity_matrix<T> >
  : DenseMatrix
  {};

} // namespace glas2


#endif
