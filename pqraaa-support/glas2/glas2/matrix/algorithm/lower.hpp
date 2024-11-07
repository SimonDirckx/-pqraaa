//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_lower_hpp
#define glas2_matrix_algorithm_lower_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename M>
  class lower_view {
    public:
      typedef typename M::value_type  value_type ;
      typedef typename M::size_type   size_type ;
      typedef lower_triangular_matrix structure ;

    public:
      lower_view( M m )
      : m_( m )
      {}

    public:
      size_type num_rows() const { return m_.num_rows() ; }
      size_type num_columns() const { return m_.num_columns() ; }

      value_type& operator()( size_type i, size_type j ) {
        assert( j<=i ) ;
        return m_(i,j) ;
      }

      value_type const& operator()( size_type i, size_type j ) const {
        assert( j<=i ) ;
        return m_(i,j) ;
      }

    public:
      M matrix() const { return m_ ; }

    private:
      M m_ ;
  } ;


  template <typename M>
  struct glas_concept< lower_view<M> >
  : DenseMatrix
  {} ;

  template <typename M>
  lower_view<M> lower( M m ) {
    return lower_view<M>( m ) ;
  }

} // namespace glas2

#endif
