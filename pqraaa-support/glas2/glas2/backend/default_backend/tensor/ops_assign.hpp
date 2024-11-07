#ifndef glas2_backend_default_backend_tensor_ops_assign_hpp
#define glas2_backend_default_backend_tensor_ops_assign_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <boost/mpl/and.hpp>
#include <cassert>

namespace glas2 {

  template <typename X, typename E>
  typename std::enable_if< boost::mpl::and_< is<DenseTensor,X>, is<DenseTensor,E> >::value, X& >::type plus_assign( default_backend, X& x, E const& e ) {
/*    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    for (typename X::size_type i=0; i<m.num_rows(); ++i) {
      for (typename X::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) += e(i,j) ;
      }
    }*/
    return x ;
  }

  template <typename X, typename E>
  typename std::enable_if< boost::mpl::and_< is<DenseTensor,X>, is<DenseTensor,E> >::value, X& >::type minus_assign( default_backend, X& x, E const& e ) {
/*    assert( m.num_rows() == e.num_rows() ) ;
    assert( m.num_columns() == e.num_columns() ) ;
    for (typename X::size_type i=0; i<m.num_rows(); ++i) {
      for (typename X::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) -= e(i,j) ;
      }
    }*/
    return x ;
  }

  template <typename X, typename E>
  typename std::enable_if< boost::mpl::and_< is<DenseTensor,X>, is<Scalar,E> >::value, X >::type multiplies_assign( default_backend, X& x, E const& e ) {
/*    for (typename X::size_type i=0; i<m.num_rows(); ++i) {
      for (typename X::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) *= e ;
      }
    }*/
    return x ;
  }

  template <typename X, typename E>
  typename std::enable_if< boost::mpl::and_< is<DenseTensor,X>, is<Scalar,E> >::value, X >::type divides_assign( default_backend, X& x, E const& e ) {
/*    for (typename X::size_type i=0; i<m.num_rows(); ++i) {
      for (typename X::size_type j=0; j<m.num_columns(); ++j) {
        m(i,j) /= e ;
      }
    }*/
    return x ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<DenseTensor,From>::value, To& >::type plus_assign( default_backend, To& to, From const& from ) {
    assert( from.order() == 1 );
    assert( from.shape()(0) == to.size() );
    typedef static_vector<int,1> shape_type ;

    for (int i=0; i<to.size(); ++i) {
      to(i) += from(shape_type({i})) ;
    }

    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<DenseTensor,From>::value, To& >::type plus_assign( default_backend, To& to, From const& from ) {
    assert( from.order() == 2 );
    assert( from.shape()(0) == to.num_rows() );
    assert( from.shape()(1) == to.num_columns() );
    typedef static_vector<int,2> shape_type ;

    for (int i=0; i<to.num_columns(); ++i) {
      for (int j=0; j<to.num_rows(); ++j) {
        to(j,i) += from(shape_type({j,i})) ;
      }
    }

    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseVector,To>::value && is<DenseTensor,From>::value, To& >::type minus_assign( default_backend, To& to, From const& from ) {
    assert( from.order() == 1 );
    assert( from.shape()(0) == to.size() );
    typedef static_vector<int,1> shape_type ;

    for (int i=0; i<to.size(); ++i) {
      to(i) -= from(shape_type({i})) ;
    }

    return to ;
  }

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<DenseTensor,From>::value, To& >::type minus_assign( default_backend, To& to, From const& from ) {
    assert( from.order() == 2 );
    assert( from.shape()(0) == to.num_rows() );
    assert( from.shape()(1) == to.num_columns() );
    typedef static_vector<int,2> shape_type ;

    for (int i=0; i<to.num_columns(); ++i) {
      for (int j=0; j<to.num_rows(); ++j) {
        to(j,i) -= from(shape_type({j,i})) ;
      }
    }

    return to ;
  }

} // namespace glas2

#endif
