#ifndef glas2_tensor_expression_tensor_vector_multiply_hpp
#define glas2_tensor_expression_tensor_vector_multiply_hpp

#include <glas2/matrix/expression/matrix_vector_multiply.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
//#include <glas2/backend/default_backend/matrix/matrix_vector_multiply.hpp>
#include <cassert>

namespace glas2 {

  template <typename T, typename M>
  class tensor_vector_multiply {
    public:
      typedef decltype( typename T::value_type() * typename M::value_type() ) value_type ;
      typedef decltype( typename T::size_type() + typename M::size_type() )   size_type ;
      typedef typename T::shape_type                                          shape_type ;

    public:
      tensor_vector_multiply( T const& t, M const& m, int mode )
      : t_( t )
      , m_( m )
      {
        assert( t.shape()(mode)==m.num_rows() ) ;
      }

    public:
      shape_type shape() { shape_type s(t_.order()-1); s(mode_) = m_.num_columns() ; return s ; } ;

      value_type operator() ( shape_type index ) const { return inner_prod( glas2::fiber( t_, index, mode_ ), m_( glas2::all(), index(mode_)) ) ; }

    private:
      T   t_ ;
      M   m_ ;
      int mode_ ;
  } ;

  template <typename T, typename M>
  struct concept< tensor_vector_multiply<T,M> >
  : DenseTensor
  {} ;

} // namespace glas2

#endif
