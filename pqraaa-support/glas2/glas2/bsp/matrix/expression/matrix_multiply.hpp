#ifndef glas2_bsp_matrix_expression_matrix_multiply_hpp
#define glas2_bsp_matrix_expression_matrix_multiply_hpp

#include <glas2/bsp/matrix/concept/bsp_internal_dense_matrix.hpp>
#include <cassert>

namespace glas2 { namespace bsp {

  template <typename M1, typename M2>
  class matrix_multiply {
    public:
      typedef decltype( typename M1::value_type() * typename M2::value_type() ) value_type ;
      typedef decltype( typename M1::size_type() * typename M2::size_type() )   size_type ;

    public:
      matrix_multiply( M1 const& m1, M2 const& m2 )
      : m1_( m1 )
      , m2_( m2 )
      {
        assert( m1.distribution()==m2.distribution() ) ;
      }

    public:
      size_type num_rows() const { return m1_.num_rows() ; }
      size_type num_columns() const { return m2_.num_columns() ; }

      M1 const& matrix1() const { return m1_ ; }
      M2 const& matrix2() const { return m2_ ; }

    private:
      M1 m1_ ;
      M2 m2_ ;
  } ;

} } // namespace glas2::bsp

namespace glas2 {

  template <typename M1, typename M2>
  struct concept< bsp::matrix_multiply<M1,M2> >
  : bsp::BSPInternalDenseMatrix
  {} ;

} // namespace glas2

#include <glas2/bsp/backend/default/matrix/matrix_multiply.hpp>

#endif
