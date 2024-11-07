//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_transformation_hpp
#define glas2_vector_type_transformation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/type/transformation.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename V, typename Op>
  class transformation< V, Op
                        , typename std::enable_if< is<DenseVector,V>::value >::type
                        >
  {
    public:
      explicit transformation( V const& v )
      : vector_( v )
      {}

    public:
      typedef typename V::size_type                         size_type ;
      typedef decltype( Op() ( typename V::value_type() ) ) value_type ;

      size_type size() const {
        return vector_.size() ;
      }

      value_type operator() ( size_type i ) const { return op_( vector_(i) ) ; }
      auto operator() ( size_type i) { return op_( vector_(i) ) ; } // ?? REFERENCE
      
    public:
      transformation& operator=( transformation const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

      template <typename E>
      transformation& operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

    public:
      V vector() const { return vector_; }

    private:
      V  vector_ ;
      Op op_ ;
  } ;

  template <typename V, typename Op>
  struct glas_concept< transformation<V,Op>
                , typename std::enable_if< is<DenseVector,V>::value>::type
                >
  : DenseVector
  {} ;

} // namespace glas2

#endif
