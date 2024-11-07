#ifndef glas2_matrix_expression_matrix_vector_multiply_hpp
#define glas2_matrix_expression_matrix_vector_multiply_hpp

#include <glas2/type/pass_reference.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
//#include <glas2/backend/default_backend/matrix/matrix_vector_multiply.hpp>
#include <cassert>

namespace glas2 {

  template <typename M, typename V>
  class matrix_vector_multiply {
    public:
      typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;
      typedef typename pass_reference<V>::type v_ref ;
      typedef typename pass_reference<M>::type m_ref ;

    public:
      matrix_vector_multiply( M const& m, V const& v )
      : m_( m )
      , v_( v )
      {
        assert( v.size()==m.num_columns() ) ;
      }

    public:
      size_type size() const { return m_.num_rows() ; }

      value_type operator() ( size_type i ) const { return inner_prod( m_(i, glas2::all()), v_ ) ; }

      template <typename R>
      typename std::enable_if< !std::is_arithmetic<R>::value
                             , matrix_vector_multiply< typename matrix_selection< M, R, all >::result_type, V >
                             >::type operator() ( R const& r ) const {
        return matrix_vector_multiply< typename matrix_selection< M, R, all >::result_type, V >( m_(r, glas2::all()), v_ ) ;
      }

    public:
      m_ref const& matrix() const { return m_ ; }
      v_ref const& vector() const { return v_ ; }

    private:
      m_ref m_ ;
      v_ref v_ ;
  } ;

  template <typename M, typename V>
  struct glas_concept< matrix_vector_multiply<M,V> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
