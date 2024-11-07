#ifndef glas2_matrix_view_matrix_selection_hpp
#define glas2_matrix_view_matrix_selection_hpp

#include <glas2/vector/view/vector_selection.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename X>
  class vec_view {
    public:
      typedef typename X::value_type value_type ;
      typedef typename X::size_type  size_type ;

    public:
      vec_view( X x )
      : matrix_( x )
      {}

    public:
      size_type size() const { return matrix_.num_rows()*matrix_.num_columns() ; }

      value_type operator() (size_type i) const { return matrix_( selection_(i) ) ; }
      value_type& operator() (size_type i) { return matrix_( selection_(i) ) ; }

      value_type operator[] (size_type i) const { return matrix_[ selection_[i] ] ; }
      value_type& operator[] (size_type i) { return matrix_[ selection_[i] ] ; }

      template <typename S>
      vector_selection< vec_view<X>, S> operator[] ( S const& s2 ) {
        return vector_selection< vec_view<X>, S>( *this, s ) ;
      }

      template <typename S2>
      typename std::enable_if< is<DenseMatrix,S2>::value, matrix_selection< matrix_selection, S2 > >::type operator[] ( S2 const& s2 ) const {
        return matrix_selection< matrix_selection<V const,S>, S2>( *this, s2 ) ;
      }

      vec_view& operator=( vec_view const& that ) {
        assert( that.size()==size() ) ;
        matrix_ = that.matrix_ ;
        return *this ;
      }

      template <typename E>
      vec_view& operator=( E const& that ) {
        return glas2::assign( *this, that ) ;
      }

    private:
      X matrix_ ;
  } ;

  template <typename X>
  class glas_concept< matrix_selection<V,S>, typename is<DenseVector,X>::type >
  : DenseVector
  {} ;

} // namespace glas2

#endif
