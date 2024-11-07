//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/view/vector_selection.hpp>
#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/external.hpp>
#include <glas2/concept/is.hpp>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename V>
  class std_vector_wrap {
    public:
      typedef V                             std_type ;
      typedef typename std_type::value_type value_type ;
      typedef typename std_type::size_type  size_type ;

      typedef ContiguousDenseVector         concept ;

    public:
        std_vector_wrap( std_type& v )
        : v_(v)
        {}

        // Copy reference !!
        std_vector_wrap( std_vector_wrap const& that )
        : v_(that.v_)
        {}


    public:
      size_type size() const { return v_.size() ; }

      value_type const& operator[] ( size_type i ) const { return v_[i] ; }
      value_type& operator[] ( size_type i ) { return v_[i] ; }

      value_type const& operator() ( size_type i ) const { return v_[i] ; }
      value_type& operator() ( size_type i ) { return v_[i] ; }

      size_type stride() const { return 1 ; }
      value_type* storage_ptr() { return &v_.start() ; }
      value_type const* storage_ptr() const { return &v_.start() ; }

      std_vector_wrap<V>& operator=( std_vector_wrap const& that ) {
        assert( that.size()==size() ) ;
        std::copy( that.storage_ptr(), that.storage_ptr()+size(), storage_ptr() ) ;
        return *this ;
      }

      template <typename E>
      std_vector_wrap<V> operator=( E const& that ) {
        return assign( *this, that ) ;
      }

      template <typename S>
      typename std::enable_if< is<DenseVector,S>::value, vector_selection< std_vector_wrap<V>, S > >::type operator[]( S const& s ) {
        return vector_selection< std_vector_wrap<V>, S >( *this, s ) ;
      }

    private:
      std_type& v_ ;
  } ;


  template <typename V>
  struct glas_concept< std_vector_wrap<V> >
  : ContiguousDenseVector
  {};

  template <typename T, typename A>
  struct glas_concept< std::vector<T,A> >
  : External
  {};

  template <typename T, typename A>
  std_vector_wrap< std::vector<T,A> > wrap( std::vector<T,A>& v) { return std_vector_wrap< std::vector<T,A> >( v ) ; }

} // namespace glas2
