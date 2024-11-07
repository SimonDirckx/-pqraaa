#ifndef glas2_matrix_expression_matrix_multiply_hpp
#define glas2_matrix_expression_matrix_multiply_hpp

#include <glas2/matrix/expression/matrix_vector_multiply.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
//#include <glas2/backend/default_backend/matrix/matrix_vector_multiply.hpp>
#include <cassert>

namespace glas2 {

  template <typename M1, typename M2>
  class matrix_multiply {
    public:
      typedef decltype( typename M1::value_type() * typename M2::value_type() ) value_type ;
      typedef decltype( typename M1::size_type() + typename M2::size_type() )   size_type ;
      typedef any_orientation orientation ;

    public:
      matrix_multiply( M1 const& m1, M2 const& m2 )
      : m1_( m1 )
      , m2_( m2 )
      {
        assert( m1.num_columns()==m2.num_rows() ) ;
      }

    public:
      size_type num_rows() const { return m1_.num_rows() ; }
      size_type num_columns() const { return m2_.num_columns() ; }

      value_type operator() ( size_type i, size_type j ) const { return inner_prod( m1_(i, glas2::all()), m2_( glas2::all(), j) ) ; }

      M1 const& matrix1() const {return m1_ ;}
      M2 const& matrix2() const {return m2_ ;}

    private:
      M1 m1_ ;
      M2 m2_ ;
  } ;

  template <typename M1, typename M2>
  struct glas_concept< matrix_multiply<M1,M2> >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
