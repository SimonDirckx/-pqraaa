#ifndef glas2_vector_view_vector_selection_hpp
#define glas2_vector_view_vector_selection_hpp

#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename S>
  class vector_selection {
    public:
      typedef typename V::value_type value_type ;
      typedef typename S::size_type  size_type ;

    public:
      vector_selection( V v, S const& s )
      : vector_( v )
      , selection_( s )
      {}

    public:
      size_type size() const { return selection_.size() ; }

      value_type operator() (size_type i) const { return vector_( selection_(i) ) ; }
      value_type& operator() (size_type i) { return vector_( selection_(i) ) ; }

      value_type operator[] (size_type i) const { return vector_[ selection_[i] ] ; }
      value_type& operator[] (size_type i) { return vector_[ selection_[i] ] ; }

      template <typename S2>
      vector_selection< vector_selection<V,S>, S2> operator[] ( S2 const& s2 ) {
        return vector_selection< vector_selection<V,S>, S2>( *this, s2 ) ;
      }

      template <typename S2>
      typename std::enable_if< is<DenseVector,S2>::value, vector_selection< vector_selection, S2 > >::type operator[] ( S2 const& s2 ) const {
        return vector_selection< vector_selection<V const,S>, S2>( *this, s2 ) ;
      }

      vector_selection& operator=( vector_selection const& that ) {
        assert( that.size()==size() ) ;
        for (size_type i=0; i<size(); ++i) {
          (*this)[i] = that[i] ;
        }
        return *this ;
      }

      template <typename E>
      vector_selection& operator=( E const& that ) {
        return glas2::assign( *this, that ) ;
      }

    private:
      V vector_ ;
      S selection_ ;
  } ;

  template <typename V, typename S>
  class glas_concept< vector_selection<V,S> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
