//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_ttt_hpp
#define glas3_array_dense_array_algorithm_ttt_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/type/ttt_operation.hpp>

#include <initializer_list>
#include <type_traits>

namespace glas3 {

template < typename X, typename Y, typename Dx, typename Dy >
typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value && is< Array, Dx >::value && is< Array, Dy >::value, ttt_operation< X, Y, Dx, Dy > >::type
ttt ( X const& x, Y const& y, Dx const& inner_dims_x, Dy const& inner_dims_y ) { // remark: ttt of empty arrays not possible
	return ttt_operation< X, Y, Dx, Dy >( x, y, inner_dims_x, inner_dims_y ) ;
}

template < typename X, typename Y, typename Dx, typename Vy >
auto ttt ( X const& x, Y const& y, Dx const& inner_dims_x, std::initializer_list<Vy> const& inner_dims_y )
-> typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value && is< Array, Dx >::value, decltype(  ttt( x, y, inner_dims_x, dense_vector<Vy>( inner_dims_y ) ) ) >::type {
	return ttt( x, y, inner_dims_x, dense_vector<Vy>( inner_dims_y ) ) ;
}

template < typename X, typename Y, typename Vx, typename Dy >
auto ttt ( X const& x, Y const& y, std::initializer_list<Vx> const& inner_dims_x, Dy const& inner_dims_y )
-> typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value && is< Array, Dy >::value, decltype( ttt( x, y, dense_vector<Vx>( inner_dims_x ), inner_dims_y ) ) >::type {
	return ttt( x, y, dense_vector<Vx>( inner_dims_x ), inner_dims_y ) ;
}

template < typename X, typename Y, typename Vx, typename Vy >
auto ttt ( X const& x, Y const& y, std::initializer_list<Vx> const& inner_dims_x, std::initializer_list<Vy> const& inner_dims_y )
-> typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value, decltype( ttt( x, y, dense_vector<Vx>( inner_dims_x ), dense_vector<Vy>( inner_dims_y ) ) ) >::type {
	return ttt( x, y, dense_vector<Vx>( inner_dims_x ), dense_vector<Vy>( inner_dims_y ) ) ;
}

template < typename X, typename Y >
auto ttt ( X const& x, Y const& y )
-> typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value, decltype( ttt( x, y, dense_vector<typename X::ndims_type>(), dense_vector<typename Y::ndims_type>() ) ) >::type {
	return ttt( x, y, dense_vector<typename X::ndims_type>(), dense_vector<typename Y::ndims_type>() ) ;
}

} // namespace glas3

#endif
