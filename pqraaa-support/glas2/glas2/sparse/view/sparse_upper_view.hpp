//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_view_sparse_upper_view_hpp
#define glas2_sparse_view_sparse_upper_view_hpp

#include <glas2/sparse/concept/sparse_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <type_traits>
#include <cmath>

namespace glas2 {

  template <typename M>
  class sparse_upper_view {
    public:
      typedef typename M::value_type  value_type ;
      typedef typename M::size_type   size_type ;
      typedef typename M::index_type  index_type ;
      typedef upper_triangular_matrix structure ;

    public:
      sparse_upper_view( M m )
      : m_( m )
      {}

    public:
      size_type num_rows() const { return m_.num_rows() ; }
      size_type num_columns() const { return m_.num_columns() ; }

    public:
      M matrix() const { return m_ ; }

    private:
      M m_ ;
  } ;


  template <typename M>
  struct glas_concept< sparse_upper_view<M> >
  : SparseMatrix
  {} ;

} // namespace glas2

#endif
