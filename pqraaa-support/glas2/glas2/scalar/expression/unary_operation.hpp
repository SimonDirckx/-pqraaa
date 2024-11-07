//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_scalar_expression_unary_operation_hpp
#define glas2_scalar_expression_unary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/expression/unary_operation.hpp>
#include <glas2/iterative/krylov/options.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <type_traits>

namespace glas2 {

  template <typename T, typename Op>
  class unary_operation< T, Op
                        , typename std::enable_if< is<Scalar,T>::value >::type
                        >
  {
    public:
      unary_operation( T const& v )
      : scalar_( v )
      {}

    public:
      typedef decltype( Op() ( T() ) ) value_type ;

      operator value_type() const { return Op()( scalar_ ) ; }

    private:
      T  scalar_ ;
  } ;

  template <typename V, typename Op>
  struct glas_concept< unary_operation<V,Op>
                , typename std::enable_if< is<Scalar,V>::value>::type
                >
  : Scalar
  {} ;

} // namespace glas2

#endif
