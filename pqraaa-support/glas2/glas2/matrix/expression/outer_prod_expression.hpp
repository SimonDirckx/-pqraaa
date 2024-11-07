#ifndef glas2_matrix_expression_outer_prod_expression_hpp
#define glas2_matrix_expression_outer_prod_expression_hpp

//#include <glas2/backend/default_backend/matrix/outer_prod.hpp>
#include <glas2/type/pass_reference.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <cassert>

namespace glas2 {

  template <typename V1, typename V2>
  class outer_prod_expression {
    public:
      typedef decltype( typename V1::value_type() * typename V2::value_type() ) value_type ;
      typedef decltype( typename V1::size_type() + typename V2::size_type() )   size_type ;
      typedef any_orientation                                                   orientation ;

    private:
      typedef typename pass_reference<V1>::type v1_ref ;
      typedef typename pass_reference<V2>::type v2_ref ;

    public:
      outer_prod_expression( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

    public:
      size_type num_rows() const { return v1_.size() ; }
      size_type num_columns() const { return v2_.size() ; }

      value_type operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows() ) ;
        assert( j>=0 && j<num_columns() ) ;
        return v1_(i) * v2_(j) ;
      }

    public:
      v1_ref vector1() const { return v1_ ; }
      v2_ref vector2() const { return v2_ ; }

    private:
      v1_ref v1_ ;
      v2_ref v2_ ;
  } ;

  template <typename V1, typename V2>
  struct glas_concept< outer_prod_expression<V1,V2> >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
