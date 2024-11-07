//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_constant_matrix_hpp
#define glas2_matrix_type_constant_matrix_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename T>
  class constant_matrix {
    public:
      constant_matrix( I const& n, I const& m, T const& v )
      : num_rows_( n )
      , num_columns_( m )
      , value_( v )
      {}

    public:
      typedef I size_type ;
      typedef T value_type ;

      size_type num_rows() const {
        return num_rows_ ;
      }

      size_type num_columns() const {
        return num_columns_ ;
      }

      value_type operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return value_ ;
      }

    private:
      size_type  num_rows_ ;
      size_type  num_columns_ ;
      value_type value_ ;
  } ;

  template <typename I, typename T>
  struct glas_concept< constant_matrix<I,T> >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
