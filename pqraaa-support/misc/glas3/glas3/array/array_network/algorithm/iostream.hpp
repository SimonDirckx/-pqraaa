//  (C) Copyright sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_iostream_hpp
#define glas3_array_array_network_algorithm_iostream_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>
#include <glas3/array/array_network/concept/CP_tensor.hpp>
#include <glas3/array/array_network/concept/tensor_train.hpp>
#include <glas3/array/array_network/concept/Tucker_tensor.hpp>

#include <iostream>
#include <type_traits>
#include <typeinfo>

namespace glas3 {

template <typename T>
typename std::enable_if< is<TensorNetwork, T>::value, std::ostream& >::type operator<< ( std::ostream& s, T const& t ) {

	typename T::ndims_type i, j ;

	s << "concept: " ;

	if ( is<TensorNetwork, T>::value ) {
		if ( is<CPTensor, T>::value ) { s << "CPTensor" << '\n' ; }
		else if ( is<TuckerTensor, T>::value ) { s << "TuckerTensor" << '\n' ; }
		else if ( is<TensorTrain, T>::value ) { s << "TensorTrain" << '\n' ; }
		else { s << "TensorNetwork" << '\n' ; }
	}

	s << "shape: { " ;
	if ( t.shape().size() > 0 ) s << t.shape()[0] ;
	for ( i = 1; i < t.shape().size(); ++i) { s << ", " << t.shape()[i] ; }
	s << " }\n" ;

	s << "rank: { " ;
	if ( t.rank().size() > 0 ) s << t.rank()[0] ;
	for ( i = 1; i < t.rank().size(); ++i) { s << ", " << t.rank()[i] ; }
	s << " }\n" ;

	s << "size: " << t.size() << '\n' ;
	s << "ndof: " << t.ndof() << '\n' ;

	s << '\n' ;

	s << "nodes2edges:" << '\n' ;
	for ( auto nodes2edges_it: t.nodes2edges() ) {
		s << "node " << nodes2edges_it.first << ": " ;
		i = 0 ;
		for ( auto e: nodes2edges_it.second ) {
			if ( i > 0 ) s << ", " ;
			s << e ;
			++i ;
		}
		s << '\n' ;
	}

	s << '\n' ;

//	s << "nodes:" << '\n' ;
//	for ( auto nodes2arrays_it: t.nodes2arrays() ) {
//		s << "node " << nodes2arrays_it.first << '\n' ;
//		s << *nodes2arrays_it.second ;
//	}

	s << "shape nodes:" << '\n' ;
	for ( auto nodes2arrays_it: t.nodes2arrays() ) {
		s << "node " << nodes2arrays_it.first << ": { " ;
		if ( nodes2arrays_it.second->shape().size() > 0 ) s << nodes2arrays_it.second->shape()[0] ;
		for ( i = 1; i < nodes2arrays_it.second->shape().size(); ++i ) { s << ", " << nodes2arrays_it.second->shape()[i] ; }
		s << " }\n" ;
	}

	s << '\n' ;

	s << "edges2nodes:" << '\n' ;
	for ( auto edges2nodes_it: t.edges2nodes() ) {
		s << "edge " << edges2nodes_it.first << ": " ;
		i = 0 ;
		for ( auto e: edges2nodes_it.second ) {
			if ( i > 0 ) s << ", " ;
			s << "( " << std::get<0>( e ) << ", " << std::get<1>( e ) << " )" ;
			++i ;
		}
		s << '\n' ;
	}

	s << '\n' ;

	s << "edges2sizes: " << '\n' ;
	for ( auto e: t.edges2sizes() )
		s << "edge " << e.first << ": " << e.second << '\n' ;

	s << '\n' ;

	s << "inner_edges: " ;
	i = 0 ;
	for ( auto e: t.inner_edges() ) {
		if ( i > 0 ) s << ", " ;
		s << e ;
		++i ;
	}
	s << '\n' ;

	s << '\n' ;

	s << "outer_edges: " ;
	i = 0 ;
	for ( auto e: t.outer_edges() ) {
		if ( i > 0 ) s << ", " ;
		s << e ;
		++i ;
	}
	s << '\n' ;

	s << '\n' ;

	s << "contraction sequence:" << '\n' ;
	j = 0 ;
	for ( auto data: t.contraction_sequence() ) {
		s << "node contraction " << j << '\n' ;
		s << "  nodes2contract: " ;
		i = 0 ;
		for ( auto node: data.nodes2contract ) {
			if ( i > 0 ) s << ", " ;
			s << node ;
			++i ;
		}
		s << '\n' ;
		s << "  edges2contract: " ;
		i = 0 ;
		for ( auto edge: data.edges2contract ) {
			if ( i > 0 ) s << ", " ;
			s << edge ;
			++i ;
		}
		s << '\n' ;
		s << "  new_node: " << data.new_node << '\n' ;
		++j ;
	}

	s << '\n' ;

	s << "contraction sequence entry:" << '\n' ;
	j = 0 ;
	for ( auto data: t.contraction_sequence_entry() ) {
		s << "node contraction " << j << '\n' ;
		s << "  nodes2contract: " ;
		i = 0 ;
		for ( auto node: data.nodes2contract ) {
			if ( i > 0 ) s << ", " ;
			s << node ;
			++i ;
		}
		s << '\n' ;
		s << "  edges2contract: " ;
		i = 0 ;
		for ( auto edge: data.edges2contract ) {
			if ( i > 0 ) s << ", " ;
			s << edge ;
			++i ;
		}
		s << '\n' ;
		s << "  new_node: " << data.new_node << '\n' ;
		++j ;
	}

	return s ;

}

} // namespace glas3

#endif
