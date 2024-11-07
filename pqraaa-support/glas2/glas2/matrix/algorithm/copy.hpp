//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_copy_hpp
#define glas2_matrix_algorithm_copy_hpp

#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/expression/copy_expession.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename M>
  class copy_expression {
    public:
      typedef typename M::value_type  value_type ;
      typedef typename M::size_type   size_type ;
      typedef typename M::structure   structure ;

    public:
      copy_expression( M const& m )
      : m_( m )
      {}

    public:
      size_type num_rows() const { return m_.num_rows() ; }
      size_type num_columns() const { return m_.num_columns() ; }

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
  struct glas_concept< copy_expression<M> >
  : DenseMatrix
  {} ;

  template <typename M>
  copy_expression<M> copy( M const& m ) {
    return copy_expression<M>( m ) ;
  }

} // namespace glas2

#endif
