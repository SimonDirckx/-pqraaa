#ifndef glas2_tensor_expression_tensor_matrix_multiply_hpp
#define glas2_tensor_expression_tensor_matrix_multiply_hpp

#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/tensor/algorithm/fiber.hpp>
#include <cassert>

namespace glas2 {

  // Replace the size of mode 'mode' by the number of rows of the matrix.
  // The tensor is multiplied on the 'mode' side by the matrix.
  template <typename T, typename M>
  class tensor_matrix_multiply {
    public:
      typedef decltype( typename T::value_type() * typename M::value_type() ) value_type ;
      typedef decltype( typename T::size_type() + typename M::size_type() )   size_type ;
      typedef typename T::shape_type                                          shape_type ;

    public:
      tensor_matrix_multiply( T const& t, M const& m, size_type mode )
      : t_( t )
      , m_( m )
      , mode_( mode )
      , index_( t.order() )
      {
        assert( t.shape()(mode)==m.num_columns() ) ;
      }

    public:
      shape_type shape() const { shape_type s(t_.order()); s=t_.shape(); s(mode_) = m_.num_rows() ; return s ; } ;

      template <typename I>
      value_type operator() ( I index ) const {
        assert( index.size()==t_.order() ) ;
        index_ = index ;
        index_(mode_) = 0 ;
        return inner_prod( glas2::fiber( t_, index_, mode_ ), m_( index(mode_), glas2::all() ) ) ;
      }

    private:
      T                                t_ ;
      M                                m_ ;
      size_type                        mode_ ;
      mutable shared_vector<size_type> index_ ;
  } ;

  template <typename T, typename M>
  struct concept< tensor_matrix_multiply<T,M> >
  : DenseTensor
  {} ;

} // namespace glas2

#endif
