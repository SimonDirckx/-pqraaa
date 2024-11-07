//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_ttt_1D_hpp
#define glas3_array_array_network_algorithm_ttt_1D_hpp

#include <glas3/concept/is.hpp>

#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/array_network/container/tensor_network.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <glas3/array/type/indirect_array.hpp>

#include <type_traits>

namespace glas3 {

template < typename X, typename Y >
typename std::enable_if< is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value
                         && std::is_same< typename X::array_type, typename Y::array_type >::value,
                         tensor_network<typename X::array_type> >::type
ttt_1D ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {

	typedef typename X::array_type   array_type ;
	typedef typename X::ndims_type   ndims_type ;

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
	{
	auto edge_x = x.outer_edges()[inner_dim_x] ;
	auto edge_y = y.outer_edges()[inner_dim_y] ;
	assert( x.edges2sizes().at( edge_x ) == y.edges2sizes().at( edge_y ) ) ;
	map_edges_x[edge_x] = i ;
	map_edges_y[edge_y] = i ;
	inner_edges.push_back( i ) ;
	++i ;
	}
	for ( j = 0; j < inner_dim_x; ++j ) {
		auto edge = x.outer_edges()[j] ;
		map_edges_x[edge] = i ;
		outer_edges.push_back( i ) ;
		++i ;
	}
	for ( j = 0; j < y.outer_edges().size(); ++j ) {
		if ( j != inner_dim_y ) {
			auto edge = y.outer_edges()[j] ;
			map_edges_y[edge] = i ;
			outer_edges.push_back( i ) ;
			++i ;
		}
	}
	for ( j = inner_dim_x + 1; j < x.outer_edges().size(); ++j ) {
		auto edge = x.outer_edges()[j] ;
		map_edges_x[edge] = i ;
		outer_edges.push_back( i ) ;
		++i ;
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

template < typename X, typename Y >
auto ttt_1D ( X const& x, Y const& y )
-> typename std::enable_if< is< TensorNetwork, X >::value && is< TensorNetwork, Y >::value, decltype ( ttt_1D( x, y, x.shape().size() - 1, 0 ) ) >::type {
	return ttt_1D( x, y, x.shape().size() - 1, 0 ) ;
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
//decltype( ttt( x, tensor_network<Y>(), inner_dims_x, inner_dims_y ) ) >::type
//{
//
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

} // namespace glas3

#endif
