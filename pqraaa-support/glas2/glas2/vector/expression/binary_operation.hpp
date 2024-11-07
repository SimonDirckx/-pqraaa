//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_expression_binary_operation_hpp
#define glas2_vector_expression_binary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/type/pass_reference.hpp>
#include <glas2/expression/binary_operation.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename V, typename Op>
  class binary_operation< S, V, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<DenseVector,V>::value>::type
                        >
  {
    public:
      binary_operation( S const& s, V const& v )
      : scalar_(s)
      , vector_( v )
      {}

      binary_operation( S const& s, V const& v, Op op )
      : scalar_(s)
      , vector_( v )
      , op_( op )
      {}

    public:
      typedef typename std::decay<V>::type                            v_type ;
      typedef typename v_type::size_type                              size_type ;
      //typedef decltype( Op() ( S(), typename v_type::value_type() ) ) value_type ;
      typedef typename std::common_type< S, typename v_type::value_type >::type value_type ;

      size_type size() const {
        return vector_.size() ;
      }
      value_type operator() ( size_type i) const { return Op()( scalar_, vector_(i) ) ; }
      //value_type operator() ( size_type i) const { return scalar_ * vector_(i) ; }

    public:
      typedef typename pass_reference<v_type>::type v_ref ;
      S const& scalar() const { return scalar_ ; }
      v_ref vector() const { return vector_ ; }

    private:
      S        scalar_ ;
      v_ref    vector_ ;
      //v_type    vector_ ;
      Op       op_ ;
  } ;

  template <typename S, typename V, typename Op>
  struct glas_concept< binary_operation<S,V,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<DenseVector,V>::value>::type
                >
  : DenseVector
  {} ;



  template <typename V, typename S, typename Op>
  class binary_operation< V, S, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<DenseVector,V>::value>::type
                        >
  {
    public:
      binary_operation( V const& v, S const& s )
      : vector_( v )
      , scalar_(s)
      {}

    public:
      typedef typename std::decay<V>::type v_type ;
      typedef typename v_type::size_type                              size_type ;
      typedef decltype( Op() ( typename v_type::value_type(), S() ) ) value_type ;

      size_type size() const {
        return vector_.size() ;
      }
      value_type operator() ( size_type i) const { return Op()( vector_(i), scalar_ ) ; }
      //auto operator() ( size_type i) const -> decltype( this->op_( this->vector_(0), this->scalar_ ) )  { return op_( vector_(i), scalar_ ) ; }

    public:
      S const& scalar() const { return scalar_ ; }
      V const& vector() const { return vector_ ; }

    private:
      typedef typename pass_reference<v_type>::type v_ref ;
      v_ref    vector_ ;
      //v_type    vector_ ;
      S        scalar_ ;
  } ;

  template <typename V, typename S, typename Op>
  struct glas_concept< binary_operation<V,S,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<DenseVector,V>::value>::type
                >
  : DenseVector
  {} ;


  template <typename V1, typename V2, typename Op>
  class binary_operation< V1, V2, Op
                        , typename std::enable_if< is<DenseVector,V1>::value && is<DenseVector,V2>::value>::type
                        >
  {
    public:
      binary_operation( V1 const& v1, V2 const& v2 )
      : v1_( v1 )
      , v2_( v2 )
      {}

    public:
      typedef typename std::decay<V1>::type v1_type ;
      typedef typename std::decay<V2>::type v2_type ;

      typedef typename v1_type::size_type                                                         size_type ;
      //typedef decltype( Op() ( typename v1_type::value_type(), typename v2_type::value_type() ) ) value_type ;
      typedef typename std::common_type< typename v1_type::value_type, typename v2_type::value_type >::type value_type ;

      size_type size() const {
        assert( v1_.size()==v2_.size() ) ;
        return v1_.size() ;
      }
      value_type operator() ( size_type i) const { return Op()( v1_(i), v2_(i) ) ; }
      value_type operator[] ( size_type i) const { return Op()( v1_[i], v2_[i] ) ; }

    private:
      typedef typename pass_reference<v1_type>::type v1_ref ;
      typedef typename pass_reference<v2_type>::type v2_ref ;
      //v1_ref  v1_ ;
      //v2_ref  v2_ ;
      v1_type  v1_ ;
      v2_type  v2_ ;
  } ;

  template <typename V1, typename V2, typename Op>
  struct glas_concept< binary_operation<V1,V2,Op>
                , typename std::enable_if< is<DenseVector,V1>::value && is<DenseVector,V2>::value >::type
                >
  : DenseVector
  {} ;


} // namespace glas2

#endif
