//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_container_tensor_network_hpp
#define glas3_array_array_network_container_tensor_network_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/type/shape_index.hpp>

#include <glas3/array/algorithm/indexing.hpp>
#include <glas3/array/array_network/algorithm/greedy_sequence.hpp>
#include <glas3/array/array_network/algorithm/full.hpp>

#include <initializer_list>
#include <vector>
#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>
#include <tuple>

#include <iostream>

namespace glas3 {

template <typename X>
class tensor_network {
	static_assert( is<Array, X>::value, "X should be an Array" ) ;

public:
	typedef X                                                   array_type ;
	typedef typename array_type::value_type                     value_type ;
	typedef typename array_type::ndims_type                     ndims_type ;
	typedef typename array_type::size_type                      size_type ;
	typedef dense_vector<size_type>                             shape_type ;
	typedef decltype( greedy_sequence( std::map<ndims_type, std::vector<ndims_type>>(), std::map<ndims_type, std::set<std::tuple<ndims_type, ndims_type>>>(),
			std::map<ndims_type, size_type>(), std::vector<ndims_type>(), std::vector<ndims_type>() ) )
			contraction_sequence_type ;
	typedef decltype( greedy_sequence( std::map<ndims_type, std::vector<ndims_type>>(), std::map<ndims_type, std::set<std::tuple<ndims_type, ndims_type>>>(),
			std::map<ndims_type, size_type>(), std::vector<ndims_type>(), std::vector<ndims_type>() ) )
			contraction_sequence_entry_type ;

public:
	tensor_network ( )
	: nodes2arrays_( boost::make_shared<std::map<ndims_type, boost::shared_ptr<array_type>>>( ) )
	, nodes2edges_( boost::make_shared<std::map<ndims_type, std::vector<ndims_type>>>( ) )
	, edges2nodes_( boost::make_shared<std::map<ndims_type, std::set<std::tuple<ndims_type, ndims_type>>>>( ) )
	, edges2sizes_( boost::make_shared<std::map<ndims_type, size_type>>( ) )
	, inner_edges_( boost::make_shared<std::vector<ndims_type>>( inner_edges ) )
	, outer_edges_( boost::make_shared<std::vector<ndims_type>>( outer_edges ) )
	, contraction_sequence_( )
    , contraction_sequence_entry_( )
	, index_( boost::make_shared<shape_index<size_type>>( ) )
	, shape_( boost::make_shared<shape_type>( ) )
	, rank_( boost::make_shared<shape_type>( ) )
	, ndof_( 0 )
	{}

	tensor_network ( std::map<ndims_type, boost::shared_ptr<array_type>> const& nodes2arrays, std::map<ndims_type, std::vector<ndims_type>> const& nodes2edges,
			std::vector<ndims_type> const& inner_edges, std::vector<ndims_type> const& outer_edges )
	: nodes2arrays_( boost::make_shared<std::map<ndims_type, boost::shared_ptr<array_type>>>( nodes2arrays ) )
	, nodes2edges_( boost::make_shared<std::map<ndims_type, std::vector<ndims_type>>>( nodes2edges ) )
	, edges2nodes_( )
	, edges2sizes_( )
	, inner_edges_( boost::make_shared<std::vector<ndims_type>>( inner_edges ) )
	, outer_edges_( boost::make_shared<std::vector<ndims_type>>( outer_edges ) )
	, contraction_sequence_( )
	, contraction_sequence_entry_( )
	, index_( )
	, shape_( )
	, rank_( )
	, ndof_( )
	{
		edges2nodes_ = boost::make_shared<std::map<ndims_type, std::set<std::tuple<ndims_type, ndims_type>>>>( ) ;
		edges2sizes_ = boost::make_shared<std::map<ndims_type, size_type>>( ) ;
		for ( auto nodes2edges_it: *nodes2edges_ ) {
			auto node = nodes2edges_it.first ;
			for ( typename array_type::ndims_type k = 0; k < nodes2edges_it.second.size(); ++k ) {
				auto edge = nodes2edges_it.second[k] ;
				(*edges2nodes_)[edge].insert( std::tuple<ndims_type, ndims_type>( node, k ) ) ;
				if ( edges2sizes_->count( edge ) ) { assert( edges2sizes_->at( edge ) == nodes2arrays_->at( node )->shape()[k] ) ; }
				else { (*edges2sizes_)[edge] = nodes2arrays_->at( node )->shape()[k] ; }
			}
		}
		contraction_sequence_ = boost::make_shared<contraction_sequence_type>( greedy_sequence( *nodes2edges_, *edges2nodes_, *edges2sizes_, *inner_edges_, *outer_edges_ ) ) ;
		contraction_sequence_entry_ = boost::make_shared<contraction_sequence_entry_type>( greedy_sequence_entry( *nodes2edges_, *edges2nodes_, *edges2sizes_, *inner_edges_, *outer_edges_ ) ) ;
		shape_ = boost::make_shared<shape_type>( no_init(), outer_edges_->size() ) ;
		for ( ndims_type k = 0; k < outer_edges_->size(); ++k ) {
			(*shape_)[k] = edges2sizes_->at( outer_edges_->at( k ) ) ;
		}
		rank_ = boost::make_shared<shape_type>( no_init(), inner_edges_->size() ) ;
		for ( ndims_type k = 0; k < inner_edges_->size(); ++k ) {
			(*rank_)[k] = edges2sizes_->at( inner_edges_->at( k ) ) ;
		}
		ndof_ = 0 ;
		for ( ndims_type k = 0; k < nodes2arrays_->size(); ++k ) {
			ndof_ += nodes2arrays_->at( k )->ndof() ;
		}
		index_ = boost::make_shared<shape_index<size_type>>( *shape_ ) ;
	}

public:
	// Copy constructor -> deep copy
	tensor_network ( tensor_network const& that )
	: nodes2arrays_( boost::make_shared<std::map<ndims_type, boost::shared_ptr<array_type>>>( ) )
	, nodes2edges_( boost::make_shared<std::map<ndims_type, std::vector<ndims_type>>>( *that.nodes2edges_ ) )
	, edges2nodes_( boost::make_shared<std::map<ndims_type, std::set<std::tuple<ndims_type, ndims_type>>>>( *that.edges2nodes_ ) )
	, edges2sizes_( boost::make_shared<std::map<ndims_type, size_type>>( *that.edges2sizes_ ) )
	, inner_edges_( boost::make_shared<std::vector<ndims_type>>( *that.inner_edges_  ) )
	, outer_edges_( boost::make_shared<std::vector<ndims_type>>( *that.outer_edges_  ) )
	, contraction_sequence_( that.contraction_sequence_ ) // -> this is not a deep copy!
    , contraction_sequence_entry_( that.contraction_sequence_entry_ )
	, index_( boost::make_shared<shape_index<size_type>>( *that.index_ ) )
	, shape_( boost::make_shared<shape_type>( *that.shape_ ) )
	, rank_( boost::make_shared<shape_type>( *that.rank_ ) )
	, ndof_( that.ndof_ )
	{}

	// Move constructor
	tensor_network ( tensor_network && ) = default ;

public:
	// Copy assignment -> do not allow
	tensor_network& operator= ( tensor_network const& that ) = delete ; // not assignable

	// Move assignment
	tensor_network& operator= ( tensor_network&& that ) = delete ;

public:
	// Constructor used in shallow_copy
	tensor_network ( boost::shared_ptr<std::map<ndims_type, boost::shared_ptr<array_type>>> nodes2arrays, boost::shared_ptr<std::map<ndims_type, std::vector<ndims_type>>> nodes2edges,
			boost::shared_ptr<std::map<ndims_type, std::set<std::tuple<ndims_type, ndims_type>>>> edges2nodes, boost::shared_ptr<std::map<ndims_type, size_type>> edges2sizes,
			boost::shared_ptr<std::vector<ndims_type>> inner_edges, boost::shared_ptr<std::vector<ndims_type>> outer_edges,
			boost::shared_ptr<contraction_sequence_type> contraction_sequence, boost::shared_ptr<contraction_sequence_entry_type> contraction_sequence_entry,
			boost::shared_ptr<shape_index<size_type>> index,
			boost::shared_ptr<shape_type> shape, boost::shared_ptr<shape_type> rank, size_type ndof)
	: nodes2arrays_( nodes2arrays )
	, nodes2edges_( nodes2edges )
	, edges2nodes_( edges2nodes )
	, edges2sizes_( edges2sizes )
	, inner_edges_( inner_edges )
	, outer_edges_( outer_edges )
    , contraction_sequence_( contraction_sequence )
	, contraction_sequence_entry_( contraction_sequence_entry )
	, index_( index )
	, shape_( shape )
	, rank_( rank )
	, ndof_( ndof )
	{}

public:
	boost::shared_ptr<std::map<ndims_type, boost::shared_ptr<array_type>>>                      nodes2arrays_ ;
	boost::shared_ptr<std::map<ndims_type, std::vector<ndims_type>>>                            nodes2edges_ ;
	boost::shared_ptr<std::map<ndims_type, std::set<std::tuple<ndims_type, ndims_type>>>>       edges2nodes_ ;
	boost::shared_ptr<std::map<ndims_type, size_type>>                                          edges2sizes_ ;
	boost::shared_ptr<std::vector<ndims_type>>                                                  inner_edges_ ;
	boost::shared_ptr<std::vector<ndims_type>>                                                  outer_edges_ ;
	boost::shared_ptr<contraction_sequence_type>                                                contraction_sequence_ ;
	boost::shared_ptr<contraction_sequence_entry_type>                                          contraction_sequence_entry_ ;
	boost::shared_ptr<shape_index<size_type>>                                                   index_ ;
	boost::shared_ptr<shape_type>                                                               shape_ ;
	boost::shared_ptr<shape_type>                                                               rank_ ;
	size_type                                                                                   ndof_ ;

public:
	tensor_network shallow_copy () const {
		return tensor_network( nodes2arrays_, nodes2edges_, edges2nodes_, edges2sizes_, inner_edges_, outer_edges_,
				contraction_sequence_, contraction_sequence_entry_, index_, shape_, rank_, ndof_ ) ;
	}

public:
	shape_type const& shape () const {
		return *shape_ ;
	}

	size_type size () const {
		size_type sz = 1 ;
		for ( ndims_type k = 0; k < shape_->size(); ++k ) {
			sz *= (*shape_)[k] ;
		}
		return sz ;
	}

	size_type ndof () const {
		return ndof_ ;
	}

public:
	shape_type rank () const {
		return *rank_ ;
	}

	std::map<typename array_type::ndims_type, boost::shared_ptr<array_type>> const& nodes2arrays () const {
		return *nodes2arrays_ ;
	}

	std::map<typename array_type::ndims_type, std::vector<typename array_type::ndims_type>>  const& nodes2edges () const {
		return *nodes2edges_ ;
	}

	std::map<typename array_type::ndims_type, std::set<std::tuple<typename array_type::ndims_type, typename array_type::ndims_type>>> const& edges2nodes () const {
		return *edges2nodes_ ;
	}

	std::map<typename array_type::ndims_type, typename array_type::size_type> const& edges2sizes () const {
		return *edges2sizes_ ;
	}

	std::vector<typename array_type::ndims_type> const& inner_edges () const {
		return *inner_edges_ ;
	}

	std::vector<typename array_type::ndims_type> const& outer_edges () const {
		return *outer_edges_ ;
	}

	std::vector<typename array_type::ndims_type>& outer_edges () {
		return *outer_edges_ ;
	}

	contraction_sequence_type const& contraction_sequence () const {
		return *contraction_sequence_ ;
	}

	contraction_sequence_entry_type const& contraction_sequence_entry () const {
		return *contraction_sequence_entry_ ;
	}

public:
	value_type	operator[] ( size_type const& i ) const {
		index_->set_lin_in( i ) ;
		return full_entry( *this, *index_ ) ;
	}

	template <typename J>
	typename std::enable_if< is<Vector, J>::value, value_type >::type operator() ( J const& j ) const {
		return full_entry( *this, j ) ;
	}

} ;

template <typename X>
struct concept< tensor_network<X>, typename std::enable_if< is<Array, X>::value >::type >
: TensorNetwork
  {};

} ; // namespace glas3

#endif
