//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_ttt_1D_hpp
#define glas3_array_dense_array_algorithm_ttt_1D_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/type/ttt_1D_operation.hpp>

#include <type_traits>

namespace glas3 {

template < typename X, typename Y >
typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value, ttt_1D_operation< X, Y > >::type
ttt_1D ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) { // remark: ttt_1D of empty arrays not possible
	return ttt_1D_operation< X, Y >( x, y, inner_dim_x, inner_dim_y ) ;
}

template < typename X, typename Y >
auto ttt_1D ( X const& x, Y const& y )
-> typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value, decltype ( ttt_1D( x, y, x.shape().size() - 1, 0 ) ) >::type {
	return ttt_1D( x, y, x.shape().size() - 1, 0 ) ;
}

} // namespace glas3

#endif
