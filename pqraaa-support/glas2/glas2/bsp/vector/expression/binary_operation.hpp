//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_vector_expression_binary_operation_hpp
#define glas2_bsp_vector_expression_binary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/expression/binary_operation.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename V, typename Op>
  class binary_operation< S, V, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPVector,V>::value>::type
                        >
  {
    public:
      typedef binary_operation< S, typename V::local_type, Op > local_type ;
      typedef typename V::distribution_type                     distribution_type ;

    public:
      binary_operation( S const& s, V const& v )
      : local_( s, v.local() )
      , distribution_( v.distribution() )
      {}

    public:
      distribution_type const& distribution() const { return distribution_ ; }

      local_type const& local() const { return local_ ; }

    private:
      local_type        local_ ;
      distribution_type distribution_ ;
  } ;

  template <typename S, typename V, typename Op>
  struct concept< binary_operation<S,V,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPVector,V>::value>::type
                >
  : bsp::BSPVector
  {} ;



  template <typename V, typename S, typename Op>
  class binary_operation< V, S, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPVector,V>::value>::type
                        >
  {
    public:
      typedef binary_operation< typename V::local_type, S, Op > local_type ;
      typedef typename V::distribution_type                     distribution_type ;

    public:
      binary_operation( V const& v, S const& s )
      : local_( v.local(), s )
      , distribution_( v.distribution() )
      {}

    public:
      distribution_type const& distribution() const { return distribution_ ; }

      local_type const& local() const { return local_ ; }

    private:
      local_type        local_ ;
      distribution_type distribution_ ;
  } ;

  template <typename V, typename S, typename Op>
  struct concept< binary_operation<V,S,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPVector,V>::value>::type
                >
  : bsp::BSPVector
  {} ;


  template <typename V1, typename V2, typename Op>
  class binary_operation< V1, V2, Op
                        , typename std::enable_if< is<bsp::BSPVector,V1>::value && is<bsp::BSPVector,V2>::value>::type
                        >
  {
    public:
      typedef binary_operation< typename V1::local_type, typename V2::local_type, Op > local_type ;
      typedef typename V1::distribution_type                                           distribution_type ;

    public:
      binary_operation( V1 const& v1, V2 const& v2 )
      : local_( v1.local(), v2.distribution() )
      , distribution_( v1.distribution() )
      {
        assert( distribution_==v2.distribution() ) ;
      }

    public:
      distribution_type const& distribution() const { return distribution_ ; }

      local_type const& local() const { return local_ ; }

    private:
      local_type        local_ ;
      distribution_type distribution_ ;
  } ;

  template <typename V1, typename V2, typename Op>
  struct concept< binary_operation<V1,V2,Op>
                , typename std::enable_if< is<bsp::BSPVector,V1>::value && is<bsp::BSPVector,V2>::value >::type
                >
  : bsp::BSPVector
  {} ;


} // namespace glas2

#endif
