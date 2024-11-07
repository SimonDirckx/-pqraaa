//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_ttt_hpp
#define glas3_array_array_network_algorithm_ttt_hpp

#include <glas3/concept/is.hpp>

#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/array_network/container/tensor_network.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <glas3/array/type/indirect_array.hpp>

#include <initializer_list>
#include <type_traits>

namespace glas3 {

template < typename X, typename Y, typename Dx, typename Dy >
typename std::enable_if< is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value
                         && std::is_same< typename X::array_type, typename Y::array_type >::value
                         && is< Array, Dx >::value && is< Array, Dy >::value, tensor_network<typename X::array_type> >::type
ttt ( X const& x, Y const& y, Dx const& inner_dims_x, Dy const& inner_dims_y ) { // remark: ttt of empty arrays not possible

	typedef typename X::array_type   array_type ;
	typedef typename X::ndims_type   ndims_type ;

	assert( inner_dims_x.size() == inner_dims_y.size() ) ;

	std::map<ndims_type, ndims_type> map_edges_x, map_edges_y ;
	std::map<ndims_type, std::vector<ndims_type>> nodes2edges ;
	std::map<ndims_type, boost::shared_ptr<array_type>> nodes2arrays ;
	std::vector<ndims_type> inner_edges, outer_edges ;
	ndims_type i, j ;

	// renumbering of old edges
	i = 0 ;
	for ( auto edge: x.inner_edges() ) {
		map_edges_x[edge] = i ;
		inner_edges.push_back( i ) ;
		++i ;
	}
	for ( auto edge: y.inner_edges() ) {
	    map_edges_y[edge] = i ;
	    inner_edges.push_back( i ) ;
	    ++i ;
	}
	for ( j = 0; j < inner_dims_x.size(); ++j ) {
		auto edge_x = x.outer_edges()[inner_dims_x[j]] ;
		auto edge_y = y.outer_edges()[inner_dims_y[j]] ;
		assert( x.edges2sizes().at( edge_x ) == y.edges2sizes().at( edge_y ) ) ;
		map_edges_x[edge_x] = i ;
		map_edges_y[edge_y] = i ;
		inner_edges.push_back( i ) ;
		++i ;
	}
	for ( auto edge: x.outer_edges() ) {
		if ( ! map_edges_x.count( edge ) ) {
			map_edges_x[edge] = i ;
			outer_edges.push_back( i ) ;
			++i ;
		}
	}
	for ( auto edge: y.outer_edges() ) {
		if ( ! map_edges_y.count( edge ) ) {
			map_edges_y[edge] = i ;
			outer_edges.push_back( i ) ;
			++i ;
		}
	}

	// construct new nodes2edges and nodes2arrays
	i = 0 ;
	for ( auto nodes2edges_it: x.nodes2edges() ) {
		auto node = nodes2edges_it.first ;
		for ( auto edge: nodes2edges_it.second ) {
			nodes2edges[i].push_back( map_edges_x.at( edge ) ) ;
		}
		nodes2arrays[i] = x.nodes2arrays().at( node ) ;
		++i ;
	}
	for ( auto nodes2edges_it: y.nodes2edges() ) {
		auto node = nodes2edges_it.first ;
		for ( auto edge: nodes2edges_it.second ) {
			nodes2edges[i].push_back( map_edges_y.at( edge ) ) ;
		}
		nodes2arrays[i] = y.nodes2arrays().at( node ) ;
		++i ;
	}

	return tensor_network<array_type>( nodes2arrays, nodes2edges, inner_edges, outer_edges ) ;
}

//template < typename X, typename Y, typename Dx, typename Dy >
//auto ttt ( X const& x, Y const& y, Dx const& inner_dims_x, Dy const& inner_dims_y )
//-> typename std::enable_if< is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value
//&& ! std::is_same< typename X::network_type::array_type, typename Y::network_type::array_type >::value
//&& is< Array, Dx >::value && is< Array, Dy >::value,
//decltype( ttt( x, tensor_network<Y>(), inner_dims_x, inner_dims_y ) ) >::type
//{}

//template < typename X, typename Y, typename Dx, typename Dy >
//auto ttt ( X const& x, Y const& y, Dx const& inner_dims_x, Dy const& inner_dims_y )
//-> typename std::enable_if< is< TensorNetwork, X >::value && is< DenseArray, Y >::value && is< Array, Dx >::value && is< Array, Dy >::value,
//decltype( ttt( x, tensor_network<typename X::array_type>, inner_dims_x, inner_dims_y ) )>::type {
//	typedef typename Y::ndims_type ndims_type ;
//
//	std::map<ndims_type, std::vector<ndims_type>> nodes2edges ;
//	std::map<ndims_type, boost::shared_ptr<Y>> nodes2arrays ;
//	std::vector<ndims_type> outer_edges( y.shape().size() ) ;
//
//	for ( ndims_type k = 0; k < y.shape().size(); ++k ) {
//		outer_edges[k] = k ;
//	}
//	nodes2edges[0] = outer_edges ;
//	nodes2arrays[0] = boost::make_shared< Y >( std::move( y.shallow_copy() ) ) ;
//
//	tensor_network<Y> y_network( nodes2arrays, nodes2edges, {}, outer_edges ) ;
//	return ttt( x, y_network, inner_dims_x, inner_dims_y ) ;
//}

//template < typename X, typename Y, typename Dx, typename Dy >
//auto ttt ( X const& x, Y const& y, Dx const& inner_dims_x, Dy const& inner_dims_y )
//-> typename std::enable_if< !is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value && is< Array, Dx >::value && is< Array, Dy >::value,
//decltype( ttt( tensor_network<X>(), y, inner_dims_x, inner_dims_y ) ) >::type
//{
//
//	typedef typename X::ndims_type ndims_type ;
//
//	std::map<ndims_type, std::vector<ndims_type>> nodes2edges ;
//	std::map<ndims_type, boost::shared_ptr<X>> nodes2arrays ;
//	std::vector<ndims_type> outer_edges( x.shape().size() ) ;
//
//	for ( ndims_type k = 0; k < x.shape().size(); ++k ) {
//		outer_edges[k] = k ;
//	}
//	nodes2edges[0] = outer_edges ;
//	nodes2arrays[0] = boost::make_shared< X >( std::move( x.shallow_copy() ) ) ;
//
//	tensor_network<X> x_network( nodes2arrays, nodes2edges, {}, outer_edges ) ;
//	return ttt( x_network, y, inner_dims_x, inner_dims_y ) ;
//}

template < typename X, typename Y, typename Dx, typename Vy >
auto ttt ( X const& x, Y const& y, Dx const& inner_dims_x, std::initializer_list<Vy> inner_dims_y )
-> typename std::enable_if< ( ( is< TensorNetwork, X >::value && is< DenseArray, Y >::value )
                           || ( is< DenseArray, X >::value && is< TensorNetwork, Y >::value )
                           || ( is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value ) )
                           && is< Array, Dx >::value, decltype(  ttt( x, y, inner_dims_x, dense_vector<Vy>( inner_dims_y ) ) ) >::type {
	return ttt( x, y, inner_dims_x, dense_vector<Vy>( inner_dims_y ) ) ;
}

template < typename X, typename Y, typename Vx, typename Dy >
auto ttt ( X const& x, Y const& y, std::initializer_list<Vx> inner_dims_x, Dy const& inner_dims_y )
-> typename std::enable_if< ( ( is< TensorNetwork, X >::value && is< DenseArray, Y >::value )
		                   || ( is< DenseArray, X >::value && is< TensorNetwork, Y >::value )
		                   || ( is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value ) )
		                   && is< Array, Dy >::value, decltype( ttt( x, y, dense_vector<Vx>( inner_dims_x ), inner_dims_y ) ) >::type {
	return ttt( x, y, dense_vector<Vx>( inner_dims_x ), inner_dims_y ) ;
}

template < typename X, typename Y, typename Vx, typename Vy >
auto ttt ( X const& x, Y const& y, std::initializer_list<Vx> inner_dims_x, std::initializer_list<Vy> inner_dims_y )
-> typename std::enable_if< ( ( is< TensorNetwork, X >::value && is< DenseArray, Y >::value )
                           || ( is< DenseArray, X >::value && is< TensorNetwork, Y >::value )
                           || ( is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value ) ),
                            decltype( ttt( x, y, dense_vector<Vx>( inner_dims_x ), dense_vector<Vy>( inner_dims_y ) ) ) >::type {
	return ttt( x, y, dense_vector<Vx>( inner_dims_x ), dense_vector<Vy>( inner_dims_y ) ) ;
}

template < typename X, typename Y >
auto ttt ( X const& x, Y const& y )
-> typename std::enable_if< ( ( is< TensorNetwork, X >::value && is< DenseArray, Y >::value )
                           || ( is< DenseArray, X >::value && is< TensorNetwork, Y >::value )
                           || ( is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value ) ),
                            decltype( ttt( x, y, dense_vector<typename X::ndims_type>(), dense_vector<typename Y::ndims_type>() ) ) >::type {
	return ttt( x, y, dense_vector<typename X::ndims_type>(), dense_vector<typename Y::ndims_type>() ) ;
}

} // namespace glas3

#endif
