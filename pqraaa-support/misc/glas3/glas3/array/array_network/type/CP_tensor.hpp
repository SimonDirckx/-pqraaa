//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_container_CP_tensor_hpp
#define glas3_array_array_network_container_CP_tensor_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/array_network/concept/CP_tensor.hpp>

#include <glas3/array/array_network/container/tensor_network.hpp>
#include <glas3/array/type/vector_wrapper.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <initializer_list>
#include <vector>
#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

#include <iostream>

namespace glas3 {

template <typename X>
class CP_tensor : public tensor_network<X> {
public:
	typedef tensor_network<X>                                  base_type ;
	typedef typename base_type::array_type                     array_type ;
	typedef typename base_type::value_type                     value_type ;
	typedef typename base_type::ndims_type                     ndims_type ;
	typedef typename base_type::size_type                      size_type ;
	typedef typename base_type::shape_type                     shape_type ;
	typedef typename base_type::contraction_sequence_type      contraction_sequence_type ;

public:
	CP_tensor ( )
	: base_type( )
	{}

	template < typename V >
	CP_tensor ( V const& value, vector_wrapper<size_type> const& shape, vector_wrapper<size_type> const& rank )
	: base_type( init_nodes2arrays( value, shape, rank ), init_nodes2edges( shape.size() ), init_inner_edges( shape.size() ), init_outer_edges( shape.size() ) )
	  {}

private:
	template < typename V, typename S, typename R >
	static std::map<ndims_type, boost::shared_ptr<X>>
	init_nodes2arrays ( V const& value, S const& shape, R const& rank ) {
		std::map<ndims_type, boost::shared_ptr<X>> nodes2arrays ;
		for ( ndims_type i = 0; i < shape.size(); ++i ) {
			nodes2arrays[i] = boost::make_shared< X >( std::move( X( value, { shape[i], rank[0] } ) ) ) ;
		}
		return nodes2arrays ;
	}

	template < typename V >
	static std::map<ndims_type, std::vector<ndims_type>>
	init_nodes2edges ( V const& ndims ) {
		std::map<ndims_type, std::vector<ndims_type>> nodes2edges ;
		for ( ndims_type i = 0; i < ndims; ++i ) {
			nodes2edges[i].push_back( i + 1 ) ;
			nodes2edges[i].push_back( 0 ) ;
		}
		return nodes2edges ;
	}

	template < typename V >
	static std::vector<ndims_type>
	init_inner_edges ( V const& ndims ) {
		std::vector<ndims_type> inner_edges ( 1, 0 ) ;
		return inner_edges ;
	}

	template < typename V >
	static std::vector<ndims_type>
	init_outer_edges ( V const& ndims ) {
		std::vector<ndims_type> outer_edges ( ndims ) ;
		for ( ndims_type i = 0; i < ndims; ++i ) {
			outer_edges[i] = i + 1 ;
		}
		return outer_edges ;
	}

public:
	// Copy constructor -> deep copy
	CP_tensor ( CP_tensor const& that )
	: base_type( that )
	{}

	// Move constructor
	CP_tensor ( CP_tensor && ) = default ;

public:
	// Copy assignment -> do not allow
	CP_tensor& operator= ( CP_tensor const& that ) = delete ; // not assignable

	// Move assignment
	CP_tensor& operator= ( CP_tensor&& that ) = delete ;

public:
	// Constructor used in shallow_copy
	CP_tensor ( base_type && base )
	: base_type( std::move( base ) )
	{}

public:
	CP_tensor shallow_copy () const {
		return CP_tensor( std::move( base_type::shallow_copy() ) ) ;
	}

} ;

template <typename X>
struct concept< CP_tensor<X>, typename std::enable_if< is<Array, X>::value >::type >
: CPTensor
  {};

} ; // namespace glas3

#endif
