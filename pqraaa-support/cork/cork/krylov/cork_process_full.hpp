//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_cork_process_full_hpp
#define cork_krylov_cork_process_full_hpp

#include <cork/options/debug_level.hpp>
#include <cork/options/backend.hpp>
#include <cork/krylov/expand_to_full.hpp>
#include <cork/krylov/cork_quadruple.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getrs.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include <cmath>
#include <type_traits>

namespace CORK { namespace krylov {

  template <typename ShiftValueType, typename Linearization, typename Quadruple>
  auto cork_process_full( Linearization& linearization, ShiftValueType const& shift, Quadruple& quadruple ) {
    typedef typename Linearization::template value_type_for<ShiftValueType> value_type ;
    glas2::matrix< value_type > A( linearization.num_rows(), linearization.num_columns() ) ;
    glas2::matrix< value_type > B( linearization.num_rows(), linearization.num_columns() ) ;

    auto fill_handle = linearization.fill_handle() ;
    fill_handle.A( A ) ;
    fill_handle.B( B ) ;

    glas2::shared_matrix< value_type > W( linearization.num_rows(), 2 ) ;
    expand_to_full( quadruple, W(glas2::all(), glas2::range(0,1)), 0 ) ;

    std::cout << "Shift = " << shift << std::endl ;
    glas2::matrix< value_type > assembled( linearization.num_rows(), linearization.num_columns() ) ;
    assembled = A - shift * B ;

    // Do one Krylov step
    auto v = W( glas2::all(), 0 ) ;
    auto w = W( glas2::all(), 1 ) ;
    w = multiply( B, v ) ;

    std::vector<int>    pivots( w.size() ) ;
    int info = boost::numeric::bindings::lapack::getrf( assembled, pivots ) ;
    assert( info>=0 ) ;
    if (info>0) throw exception::linear_solver_failure() ;

    info = boost::numeric::bindings::lapack::getrs( assembled, pivots, w ) ;
    assert( info==0 ) ;
    //std::cout << "w = " << w << std::endl ;

    glas2::shared_matrix<value_type> H(2,1) ;
    H(0,0) = inner_prod( conj(v), w ) ;
    w -= H(0,0) * v ;
    value_type g = inner_prod( conj(v), w ) ;
    w -= g * v ;
    H(0,0) += g ;
    H(1,0) = norm_2(w) ;
    w /= H(1,0) ;

    std::cout << "Ortho vectors " << inner_prod( conj(v), w ) << std::endl ;

    return std::tuple(W, H) ;
  } // cork_process_full()

} } // namespace CORK::krylov


#endif
