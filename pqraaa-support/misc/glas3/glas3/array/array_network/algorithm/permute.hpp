//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_permute_hpp
#define glas3_array_array_network_algorithm_permute_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/array_network/container/tensor_network.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/CP_tensor.hpp>

#include <type_traits>
#include <initializer_list>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

namespace glas3 {

template < typename T, typename D >
typename std::enable_if< is< TensorNetwork, T >::value && ! is< CPTensor, T >::value
                    && ! is< TensorTrain, T >::value && ! is< TuckerTensor, T >::value && is< Array, D >::value, tensor_network<typename T::array_type> >::type
permute( T const& t, D const& dim_order ) {
	assert( t.outer_edges_->size() == dim_order.size() ) ;
	auto outer_edges = boost::make_shared<std::vector<typename T::ndims_type>>( t.outer_edges_->size() ) ;
	for ( typename D::size_type i = 0; i < dim_order.size(); ++i ) {
	  	outer_edges->at(i) = t.outer_edges_->at( dim_order[i] ) ;
	}
//	auto shape = boost::make_shared<typename T::shape_type>( no_init(), outer_edges->size() ) ;
//	for ( typename T::ndims_type k = 0; k < outer_edges->size(); ++k ) {
//	    (*shape)[k] = t.edges2sizes_->at( outer_edges->at( k ) ) ;
//	}
//	auto index = boost::make_shared<shape_index<typename T::size_type>>( *shape ) ;
//   	return tensor_network<typename T::array_type>( t.nodes2arrays_, t.nodes2edges_, t.edges2nodes_, t.edges2sizes_, t.inner_edges_, outer_edges, t.contraction_sequence_, t.contraction_sequence_entry_, index, shape, t.rank_, t.ndof_ ) ;
	// Doesn't work because contraction_sequence depends on outer_edges, therefore:
   	return tensor_network<typename T::array_type>( t.nodes2arrays(), t.nodes2edges(), t.inner_edges(), *outer_edges ) ;
}

template < typename T, typename D >
typename std::enable_if< is< CPTensor, T >::value && is< Array, D >::value, CP_tensor<typename T::array_type> >::type
permute( T const& t, D const& dim_order ) {
	tensor_network<typename T::array_type> const& t_base = t ;
	return CP_tensor<typename T::array_type>( permute( t_base, dim_order ) ) ;
}

template < typename T, typename D >
typename std::enable_if< is< TensorTrain, T >::value && is< Array, D >::value, tensor_train<typename T::array_type> >::type
permute( T const& t, D const& dim_order ) {
	tensor_network<typename T::array_type> const& t_base = t ;
	return tensor_train<typename T::array_type>( permute( t_base, dim_order ) ) ;
}

template < typename T, typename D >
typename std::enable_if< is< TuckerTensor, T >::value && is< Array, D >::value, Tucker_tensor<typename T::array_type> >::type
permute( T const& t, D const& dim_order ) {
	tensor_network<typename T::array_type> const& t_base = t ;
	return Tucker_tensor<typename T::array_type>( permute( t_base, dim_order ) ) ;
}

template < typename T >
auto permute( T const& t, std::initializer_list< typename T::ndims_type > dim_order )
-> typename std::enable_if< is< TensorNetwork, T >::value, decltype( permute( t, dense_vector<typename T::ndims_type>( dim_order ) ) ) >::type {
    return permute( t, dense_vector<typename T::ndims_type>( dim_order ) ) ;
}

} // namespace glas3

#endif
