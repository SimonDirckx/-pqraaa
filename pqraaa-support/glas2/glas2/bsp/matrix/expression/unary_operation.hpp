//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_matrix_expression_unary_operation_hpp
#define glas2_bsp_matrix_expression_unary_operation_hpp

#include <glas2/bsp/matrix/concept/bsp_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/matrix/expression/unary_operation.hpp>
#include <functional>

namespace glas2 {

  template <typename V, typename Op>
  class unary_operation< V, Op
                        , typename std::enable_if< is<bsp::BSPMatrix,V>::value >::type
                        >
  {
    public:
      typedef unary_operation< typename V::local_type, Op > local_type ;

    public:
      unary_operation( V const& v )
      : local_( v.local() )
      , distribution_( v.distribution() )
      {}

      unary_operation( V const& v, Op const& )
      : local_( v.local() )
      , distribution_( v.distribution() )
      {}

      unary_operation( unary_operation const& that )
      : local_( that.local_ )
      {}

    public:
      typedef typename V::distribution_type distribution_type ;
      distribution_type const& distribution() const { return distribution_ ; }

      local_type const& local() const { return local_ ; }

    private:
      local_type        local_ ;
      distribution_type distribution_ ;
  } ;

  template <typename V, typename F>
  struct concept< glas2::unary_operation<V,F>, typename std::enable_if< is<bsp::BSPMatrix,V>::value >::type >
  : bsp::BSPMatrix
  {} ;

} // namespace glas2

#endif
