//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_type_Tucker_tensor_hpp
#define glas3_array_array_network_type_Tucker_tensor_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/array_network/concept/Tucker_tensor.hpp>

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
#include <map>

#include <iostream>

namespace glas3 {

template <typename X>
class Tucker_tensor : public tensor_network<X> {
public:
	typedef tensor_network<X>                                  base_type ;
	typedef typename base_type::array_type                     array_type ;
	typedef typename base_type::value_type                     value_type ;
	typedef typename base_type::ndims_type                     ndims_type ;
	typedef typename base_type::size_type                      size_type ;
	typedef typename base_type::shape_type                     shape_type ;
	typedef typename base_type::contraction_sequence_type      contraction_sequence_type ;

public:
	Tucker_tensor ( )
    : base_type( )
    {}

	template < typename V >
	Tucker_tensor ( V const& value, vector_wrapper<size_type> const& shape, vector_wrapper<size_type> const& rank )
	: base_type( init_nodes2arrays( value, shape, rank ), init_nodes2edges( shape.size() ), init_inner_edges( shape.size() ), init_outer_edges( shape.size() ) )
	  {}

private:
	template < typename V, typename S, typename R >
	static std::map<ndims_type, boost::shared_ptr<X>>
	init_nodes2arrays ( V const& value, S const& shape, R const& rank ) {
		std::map<ndims_type, boost::shared_ptr<X>> nodes2arrays ;
		ndims_type i = 0 ;
		nodes2arrays[i] = boost::make_shared< X >( std::move( X( value, rank ) ) ) ;
		for ( i = 0; i < shape.size(); ++i ) {
			nodes2arrays[i + 1] = boost::make_shared< X >( std::move( X( value, { shape[i], rank[i] } ) ) ) ;
		}
		return nodes2arrays ;
	}

	template < typename V >
	static std::map<ndims_type, std::vector<ndims_type>>
	init_nodes2edges ( V const& ndims ) {
		std::map<ndims_type, std::vector<ndims_type>> nodes2edges ;
		ndims_type i = 0, j ;
		for ( j = 0; j < ndims; ++j ) {
			nodes2edges[i].push_back( j ) ;
		}
		for ( i = 0; i < ndims; ++i ) {
			nodes2edges[i + 1].push_back( i + ndims ) ;
			nodes2edges[i + 1].push_back( i ) ;
		}
		return nodes2edges ;
	}

	template < typename V >
	static std::vector<ndims_type>
	init_inner_edges ( V const& ndims ) {
		std::vector<ndims_type> inner_edges ( ndims ) ;
		for ( ndims_type i = 0; i < ndims; ++i ) {
			inner_edges[i] = i ;
		}
		return inner_edges ;
	}

	template < typename V >
	static std::vector<ndims_type>
	init_outer_edges ( V const& ndims ) {
		std::vector<ndims_type> outer_edges ( ndims ) ;
		for ( ndims_type i = 0; i < ndims; ++i ) {
			outer_edges[i] = i + ndims ;
		}
		return outer_edges ;
	}

public:
	// Copy constructor -> deep copy
	Tucker_tensor ( Tucker_tensor const& that )
	: base_type( that )
	{}

	// Move constructor
	Tucker_tensor ( Tucker_tensor && ) = default ;

public:
	// Copy assignment -> do not allow
	Tucker_tensor& operator= ( Tucker_tensor const& that ) = delete ; // not assignable

	// Move assignment
	Tucker_tensor& operator= ( Tucker_tensor&& that ) = delete ;

public:
	// Constructor used in shallow_copy
	Tucker_tensor ( base_type && base )
    : base_type( std::move( base ) )
    {}

public:
	Tucker_tensor shallow_copy () const {
		return Tucker_tensor( std::move( base_type::shallow_copy() ) ) ;
	}

} ;

template <typename X>
struct concept< Tucker_tensor<X>, typename std::enable_if< is<Array, X>::value >::type >
: TuckerTensor
  {};

} ; // namespace glas3

#endif
