#ifndef glas2_bsp_matrix_expression_matrix_vector_multiply_hpp
#define glas2_bsp_matrix_expression_matrix_vector_multiply_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <cassert>
//#include <iostream>

namespace glas2 { namespace bsp {

  template <typename M, typename V>
  class matrix_vector_multiply {
    public:
      typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;

    public:
      matrix_vector_multiply( M const& m, V const& v )
      : m_( m )
      , v_( v )
      {
	//std::cout<<"coucou"<<std::endl;
        assert( v.size()==m.num_columns() ) ;
      }

    public:
      size_type size() const { return m_.num_rows() ; }

      value_type operator() ( size_type i ) const { return inner_prod( m_(i, glas2::all()), v_ ) ; }

    public:
      M const& matrix() const { return m_ ; }
      V const& vector() const { return v_ ; }

    private:
      M m_ ;
      V v_ ;
  } ;

} } // namespace glas2::bsp

namespace glas2 {

  template <typename M, typename V>
  struct concept< bsp::matrix_vector_multiply<M,V> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
