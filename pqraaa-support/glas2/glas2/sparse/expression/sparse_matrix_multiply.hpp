#ifndef glas2_sparse_expression_sparse_matrix_multiply_hpp
#define glas2_sparse_expression_sparse_matrix_multiply_hpp

#include <glas2/matrix/concept/matrix.hpp>
#include <cassert>

namespace glas2 {

  template <typename S, typename M>
  class sparse_matrix_multiply {
    public:
      typedef decltype( typename S::value_type() * typename M::value_type() ) value_type ;
      typedef decltype( typename S::size_type() + typename M::size_type() )   size_type ;

    public:
      sparse_matrix_multiply( S const& s, M const& m )
      : s_( s )
      , m_( m )
      {
        assert( s.num_columns()==m.num_rows() ) ;
      }

    public:
      size_type num_rows() const { return s_.num_rows() ; }
      size_type num_columns() const { return m_.num_columns() ; }

    public:
      S const& sparse() const { return s_ ; }
      M const& matrix() const { return m_ ; }

    private:
      S s_ ;
      M m_ ;
  } ;

  template <typename S, typename M>
  struct glas_concept< sparse_matrix_multiply<S,M> >
  : Matrix
  {} ;

} // namespace glas2


#include <glas2/algorithm/ops_assign.hpp>
#include <glas2/backend/default_backend/sparse/sparse_matrix_multiply.hpp>

namespace glas2 {

  template <typename To, typename S, typename E>
  To assign( To to, sparse_matrix_multiply<S,E> const& from ) {
    return assign( current_backend(), to, from ) ;
  }

} // namespace glas2

#endif
