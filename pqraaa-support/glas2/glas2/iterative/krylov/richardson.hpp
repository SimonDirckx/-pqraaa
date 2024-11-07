//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_iterative_krylov_richardson_hpp
#define glas_toolbox_iterative_krylov_richardson_hpp

#include <glas/config.hpp>
#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

#include <glas/toolbox/iterative/linear_operator/category.hpp>
#include <glas/toolbox/iterative/krylov/tolerance.hpp>
#include <glas/container/dense_vector.hpp>
#include <glas/algorithm/norm_2.hpp>
#include <glas/concept/value_type.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <cassert>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace glas {

  template <typename X, typename Y, typename Prec, typename Op, typename Report>
  void richardson( Op const& op, Prec const& prec, X& x, Y& y, Report& report, options const& opt ) {
    BOOST_STATIC_ASSERT( (boost::is_same< left_preconditioned_operator_tag, typename linear_operator_category<Op>::type >::value) ) ;
    typedef typename value_type<X>::type  value_type ;
    typedef dense_vector< value_type >    container_type ;

    assert( size(x) == size(y) ) ;

    container_type residu( size(x ) ) ;
    prec( y ) ;

    op( x, residu ) ;
    residu -= y ;
    double res_norm_0 = norm_2( residu ) ;
    report( res_norm_0, 0 ) ;

    for ( unsigned int it=1 ; it<opt.max_mat_vec_ ; ++it ) {
      x -= residu ;
      op( x, residu ) ;
      residu -= y ;

      double res_norm = norm_2( residu ) ;

      report( res_norm, it ) ;
      if ( res_norm < tolerance( opt, res_norm_0 ) ) return ;
    }
  } // richardson
}

#endif
