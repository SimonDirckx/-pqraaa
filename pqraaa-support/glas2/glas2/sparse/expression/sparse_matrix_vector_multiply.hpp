#ifndef glas2_sparse_expression_sparse_matrix_vector_multiply_hpp
#define glas2_sparse_expression_sparse_matrix_vector_multiply_hpp

#include <glas2/backend/current_backend.hpp>
#include <glas2/sparse/algorithm/inner_prod.hpp>
#include <glas2/sparse/concept/compressed_sparse_matrix.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename V>
  class sparse_matrix_vector_multiply {
    public:
      typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;

    public:
      sparse_matrix_vector_multiply( M const& m, V const& v )
      : m_( m )
      , v_( v )
      {
        assert( v.size()==m.num_columns() ) ;
      }

    public:
      size_type size() const { return m_.num_rows() ; }

      value_type operator()( size_type i) const {
        return inner_prod( m_( i, glas2::all() ), v_ ) ;
      }

      M const& sparse() const { return m_ ; }
      V const& vector() const { return v_ ; }

    private:
      M m_ ;
      V v_ ;
  } ;

  template <typename M, typename V>
  struct glas_concept< sparse_matrix_vector_multiply<M,V>, typename std::enable_if< is< CompressedSparseMatrix,M >::value >::type >
  : std::conditional< std::is_same< row_major, typename M::orientation >::value
                    , DenseVector
                    , Vector
                    >::type
  {} ;

  template <typename M, typename V>
  struct glas_concept< sparse_matrix_vector_multiply<M,V>, typename std::enable_if< !is< CompressedSparseMatrix,M >::value >::type >
  : Vector
  {} ;

} // namespace glas2

#include <glas2/algorithm/ops_assign.hpp>
#include <glas2/backend/default_backend/sparse/sparse_matrix_vector_multiply.hpp>

namespace glas2 {

  template <typename To, typename S, typename E>
  To assign( To to, sparse_matrix_vector_multiply<S,E> const& from ) {
    return assign( current_backend(), to, from ) ;
  }

} // namespace glas2

#endif
