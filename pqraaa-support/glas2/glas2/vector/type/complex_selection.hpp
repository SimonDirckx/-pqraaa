#ifndef glas2_vector_type_complex_selection_hpp
#define glas2_vector_type_complex_selection_hpp

#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename S>
  class complex_selection {
    public:
      typedef typename V::value_type                        complex_value_type ;
      typedef decltype( std::real( complex_value_type() ) ) value_type ;
      typedef typename V::size_type                         size_type ;

    public:
      complex_selection( V v )
      : vector_( v )
      {}

    public:
      V vector() const& { return vector_ ; }

    public:
      size_type size() const { return vector_.size() ; }

      value_type operator() (size_type i) const { return S()( vector_(i) ) ; }
      value_type operator[] (size_type i) const { return S()( vector_(i) ) ; }

      complex_selection& operator=( complex_selection const& that ) {
        glas2::assign( *this, that ) ;
        return *this ;
      }

      template <typename E>
      complex_selection& operator=( E const& that ) {
        glas2::assign( *this, that ) ;
        return *this ;
      }

    public:
      template <typename I>
      typename std::enable_if< !std::is_integral<I>::value, typename vector_selection< complex_selection, I >::result_type>::type operator()( I const& s ) {
        return vector_selection< complex_selection, I >::apply( *this, s ) ;
      }

      template <typename I>
      typename std::enable_if< !std::is_integral<I>::value, typename vector_selection< complex_selection, I >::result_type>::type operator()( I const& s ) const {
        return vector_selection< complex_selection, I >::apply( *this, s ) ;
      }

    private:
      V vector_ ;
  } ;

  template <typename V, typename S>
  struct glas_concept< complex_selection<V,S> >
  : glas_concept<V>
  {} ;

  template <typename V, typename S, typename S2>
  struct vector_selection< complex_selection<V,S>, S2 > {
    typedef typename vector_selection< V, S2 >::result_type v_result_type ;
    typedef complex_selection< v_result_type, S >           result_type ;
   
    static result_type apply( vector_selection<V,S> v, S2 r ) {
      return result_type( v.vector()(r) ) ;
    }
  } ;

} // namespace glas2

#endif
