//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_ops_hpp
#define glas3_array_array_network_algorithm_ops_hpp

#include <glas3/concept/is.hpp>

#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/blocks_array.hpp>
#include <glas3/array/type/array_wrapper.hpp>
#include <glas3/array/array_network/container/tensor_network.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/CP_tensor.hpp>

#include <type_traits>
#include <set>
#include <map>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

namespace glas3 {

template <typename A1, typename A2>
typename std::enable_if< is<TensorNetwork, A1>::value && is<TensorNetwork, A2>::value
                  && ! ( is<CPTensor, A1>::value && is<CPTensor, A2>::value )
                  && ! ( is<TensorTrain, A1>::value && is<TensorTrain, A2>::value )
                  && ! ( is<TuckerTensor, A1>::value && is<TuckerTensor, A2>::value ),
                  tensor_network<blocks_array<decltype( typename A1::value_type() + typename A2::value_type() )> > >::type
operator+( A1 const& a1, A2 const& a2 ) {
	typedef decltype( typename A1::value_type() + typename A2::value_type() )   value_type ;
	typedef blocks_array<value_type>                                            array_type ;
	typedef typename array_type::ndims_type                                     ndims_type ;

	// check whether network topologies match
	for ( auto nodes2edges_it: a1.nodes2edges() ) {
    	assert( a2.nodes2edges().count( nodes2edges_it.first ) ) ;
    	assert( a1.nodes2edges().at( nodes2edges_it.first ) == a2.nodes2edges().at( nodes2edges_it.first ) ) ;
    }
    assert( a1.outer_edges() == a2.outer_edges() ) ;

    // do summation -> cat/diag of corresponding nodes
	std::set<ndims_type> set_outer_edges( a1.outer_edges().begin(), a1.outer_edges().end() ) ;
	std::set<ndims_type> dims2diag ;
    std::map<ndims_type, boost::shared_ptr<array_type>> nodes2arrays ;
    for ( auto nodes2arrays_it: a1.nodes2arrays() ) {
    	auto node = nodes2arrays_it.first ;
    	dims2diag.clear() ;
    	std::ptrdiff_t i = 0 ;
    	for ( auto edge: a1.nodes2edges().at( node ) ) {
    		if ( ! set_outer_edges.count( edge ) ) {
    			dims2diag.insert( i ) ;
    		}
    		++i ;
    	}
    	//nodes2arrays[node] = boost::make_shared<array_type>( std::move( blocks_array<value_type>( { *a1.nodes2arrays().at( node ), *a2.nodes2arrays().at( node ) }, dims2diag ) ) ) ;
    	// DOESN'T work -> again strange bug with array_wrapper or initializer_list
    	std::vector<boost::shared_ptr<array_wrapper<value_type>>> dum ;
    	dum.push_back( boost::make_shared<array_wrapper<value_type>>( *a1.nodes2arrays().at( node ) ) ) ;
    	dum.push_back( boost::make_shared<array_wrapper<value_type>>( *a2.nodes2arrays().at( node ) ) ) ;
    	nodes2arrays[node] = boost::make_shared<array_type>( dum, dims2diag ) ;
    }

    return tensor_network<array_type>( nodes2arrays, a1.nodes2edges(), a1.inner_edges(), a1.outer_edges() ) ;
}

template <typename A1, typename A2>
typename std::enable_if< is<TensorTrain, A1>::value && is<TensorTrain, A2>::value ,
                      tensor_train<blocks_array<decltype( typename A1::value_type() + typename A2::value_type() )> > >::type
operator+( A1 const& a1, A2 const& a2 ) {
	tensor_network<typename A1::array_type> const& a1_base = a1 ;
	tensor_network<typename A2::array_type> const& a2_base = a2 ;
	return tensor_train<blocks_array<decltype( typename A1::value_type() + typename A2::value_type() )>>( a1_base + a2_base ) ;
}

template <typename A1, typename A2>
typename std::enable_if< is<CPTensor, A1>::value && is<CPTensor, A2>::value ,
                      CP_tensor<blocks_array<decltype( typename A1::value_type() + typename A2::value_type() )> > >::type
operator+( A1 const& a1, A2 const& a2 ) {
	tensor_network<typename A1::array_type> const& a1_base = a1 ;
	tensor_network<typename A2::array_type> const& a2_base = a2 ;
	return CP_tensor<blocks_array<decltype( typename A1::value_type() + typename A2::value_type() )>>( a1_base + a2_base ) ;
}

template <typename A1, typename A2>
typename std::enable_if< is<TuckerTensor, A1>::value && is<TuckerTensor, A2>::value ,
                      Tucker_tensor<blocks_array<decltype( typename A1::value_type() + typename A2::value_type() )> > >::type
operator+( A1 const& a1, A2 const& a2 ) {
	tensor_network<typename A1::array_type> const& a1_base = a1 ;
	tensor_network<typename A2::array_type> const& a2_base = a2 ;
	return Tucker_tensor<blocks_array<decltype( typename A1::value_type() + typename A2::value_type() )>>( a1_base + a2_base ) ;
}

} // namespace glas3

#endif
