//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_multiply_hpp
#define glas3_array_array_network_algorithm_multiply_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>
#include <glas3/array/array_network/concept/CP_tensor.hpp>
#include <glas3/array/array_network/concept/tensor_train.hpp>
#include <glas3/array/array_network/concept/Tucker_tensor.hpp>

#include <glas3/array/array_network/container/tensor_network.hpp>
#include <glas3/array/array_network/type/CP_tensor.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>

#include <glas2/concept/is.hpp>
#include <glas2/sparse/concept/sparse_matrix.hpp>

#include <glas3/array/dense_array/algorithm/multiply.hpp>

#include <type_traits>
#include <iostream>
#include <tuple>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <typeinfo>

namespace glas3 {

template < typename X, typename Y >
typename std::enable_if< is< TensorNetwork, X >::value && ! is< CPTensor, X >::value
                    && ! is< TensorTrain, X >::value && ! is< TuckerTensor, X >::value
                    && ( is< DenseMatrix, Y >::value || glas2::is< glas2::SparseMatrix, Y >::value )
                    && std::is_same< typename X::array_type, decltype( multiply( typename X::array_type(), Y(), 0, 0 ) ) >::value,
                    tensor_network<typename X::array_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename X::ndims_type inner_dim_y ) {
	typename X::ndims_type node, dim ;
	std::tie( node, dim ) = *x.edges2nodes().at( x.outer_edges().at( inner_dim_x ) ).begin() ;
	auto nodes2arrays = x.nodes2arrays() ;
	auto array = nodes2arrays.at( node )->shallow_copy() ;

	nodes2arrays.at( node ) = boost::make_shared<typename X::array_type>( multiply( array, y, dim, inner_dim_y ) ) ;

	if ( array.shape()[dim] == nodes2arrays.at( node )->shape()[dim] ) {
		return tensor_network<typename X::array_type>( boost::make_shared<std::map<typename X::ndims_type, boost::shared_ptr<typename X::array_type>>>( std::move( nodes2arrays ) ), x.nodes2edges_, x.edges2nodes_, x.edges2sizes_, x.inner_edges_, x.outer_edges_, x.contraction_sequence_, x.contraction_sequence_entry_, x.index_, x.shape_, x.rank_, x.ndof_ ) ;
	}
	else {
		return tensor_network<typename X::array_type>( nodes2arrays, x.nodes2edges(), x.inner_edges(), x.outer_edges() ) ;
	}
}

template < typename X, typename Y >
typename std::enable_if< is< CPTensor, X >::value
                    && ( is< DenseMatrix, Y >::value || glas2::is< glas2::SparseMatrix, Y >::value )
                    && std::is_same< typename X::array_type, decltype( multiply( typename X::array_type(), Y(), 0, 0 ) ) >::value,
                    CP_tensor<typename X::array_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename X::ndims_type inner_dim_y ) {
	tensor_network<typename X::array_type> const& t_base = x ;
    return CP_tensor<typename X::array_type>( multiply ( t_base, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< TuckerTensor, X >::value
                    && ( is< DenseMatrix, Y >::value || glas2::is< glas2::SparseMatrix, Y >::value )
                    && std::is_same< typename X::array_type, decltype( multiply( typename X::array_type(), Y(), 0, 0 ) ) >::value,
                    Tucker_tensor<typename X::array_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename X::ndims_type inner_dim_y ) {
	tensor_network<typename X::array_type> const& t_base = x ;
    return Tucker_tensor<typename X::array_type>( multiply ( t_base, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< TensorTrain, X >::value
                    && ( is< DenseMatrix, Y >::value || glas2::is< glas2::SparseMatrix, Y >::value )
                    && std::is_same< typename X::array_type, decltype( multiply( typename X::array_type(), Y(), 0, 0 ) ) >::value,
                    tensor_train<typename X::array_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename X::ndims_type inner_dim_y ) {
	tensor_network<typename X::array_type> const& t_base = x ;
    return tensor_train<typename X::array_type>( multiply ( t_base, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< TensorNetwork, X >::value
                    && ( is< DenseArray, Y >::value && ! is< DenseMatrix, Y >::value )
                    && std::is_same< typename X::array_type, decltype( multiply( typename X::array_type(), Y(), 0, 0 ) ) >::value,
                    tensor_network<typename X::array_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typename X::ndims_type node, dim, new_edge ;
	std::tie( node, dim ) = *x.edges2nodes().at( x.outer_edges().at( inner_dim_x ) ).begin() ;
	auto nodes2arrays = x.nodes2arrays() ;
	auto nodes2edges = x.nodes2edges() ;
	auto inner_edges = x.inner_edges() ;
	auto outer_edges = x.outer_edges() ;

	auto max_edge = outer_edges.at( 0 ) ;
	for ( auto edge: outer_edges ) {
		if ( edge > max_edge ) { max_edge = edge ; }
	}
	for ( auto edge: inner_edges ) {
		if ( edge > max_edge ) { max_edge = edge ; }
	}
	typename X::ndims_type n_edges2add = y.shape().size() - 1 ;
    nodes2edges.at( node ).erase( nodes2edges.at( node ).begin() + dim ) ;
    outer_edges.erase( outer_edges.begin() + inner_dim_x ) ;
    for ( typename X::ndims_type i = 0; i < n_edges2add; ++i ) {
    	nodes2edges.at( node ).insert( nodes2edges.at( node ).begin() + dim + i, max_edge + 1 + i ) ;
    	outer_edges.insert( outer_edges.begin() + inner_dim_x + i, max_edge + 1 + i ) ;
    }

	nodes2arrays.at( node ) = boost::make_shared<typename X::array_type>( multiply( *nodes2arrays.at( node ), y, dim, inner_dim_y ) ) ;

	return tensor_network<typename X::array_type>( nodes2arrays, nodes2edges, inner_edges, outer_edges ) ;
}

} // namespace glas3

#endif
