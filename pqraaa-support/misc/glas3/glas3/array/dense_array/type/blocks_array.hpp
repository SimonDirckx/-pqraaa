//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_blocks_array_hpp
#define glas3_array_type_blocks_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/array_wrapper.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>
#include <glas3/array/dense_array/container/primitive_dense_vector.hpp>
#include <glas3/array/dense_array/type/range.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/assign.hpp>
#include <glas3/array/dense_array/algorithm/block_selection.hpp>
#include <glas3/array/algorithm/indexing.hpp>

#include <initializer_list>
#include <vector>
#include <set>
#include <cassert>
#include <type_traits>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>
#include <iostream>

namespace glas3 {

template <typename V, std::ptrdiff_t ndims = -1>
class blocks_array {
	static_assert( std::is_arithmetic<V>::value || boost::is_complex<V>::value, "V should be arithmetic" ) ;

public:
	typedef V                                 value_type ;
	typedef dense_vector<std::ptrdiff_t>      shape_type ;
	typedef typename shape_type::value_type   size_type ;
	typedef typename shape_type::size_type    ndims_type ;

public:
	blocks_array ( )
	: blocks_( )
    , dims2diag_( )
	, acc_block_shapes_( )
	, shape_( boost::make_shared< shape_type >( ) )
    , dim_factor_( )
	, size_( 0 )
    , ndof_( 0 )
    {
		if ( ndims >= 0 ) {
			assert ( ndims == 0 ) ;
		}
    }

	blocks_array ( std::vector<array_wrapper<value_type>> const& blocks, std::set<ndims_type> dims2diag )
	: blocks_( boost::make_shared<std::vector<boost::shared_ptr<array_wrapper<value_type>>>>( ) )
	, dims2diag_( boost::make_shared<std::set<ndims_type>>( std::move( dims2diag ) ) )
	, acc_block_shapes_( )
	, shape_( )
	, dim_factor_( )
	, size_( 1 )
	, ndof_( 0 )
	{
		assert( dims2diag_->size() >= 1 ) ;
		for ( auto blocks_it = blocks.begin(); blocks_it != blocks.end(); ++blocks_it ){
			blocks_->push_back( boost::make_shared<array_wrapper<value_type>>( blocks_it->shallow_copy() ) ) ;
		}
		ndims_type d = blocks_->at(0)->shape().size() ;
    	if ( ndims >= 0 ) {
    		assert ( ndims == d ) ;
    	}
		ndof_ += blocks_->at(0)->size() ;
		acc_block_shapes_ = boost::make_shared<dense_matrix<size_type>>( dense_matrix<size_type>( no_init(), {ndims_type(blocks_->size()), d} ) ) ;
		(*acc_block_shapes_)({0, range<>(0, d)}) = blocks_->at(0)->shape() ;
		for ( size_type i = 1; i < blocks_->size(); ++i ) {
			assert( blocks_->at( i )->shape().size() == d ) ;
			for ( ndims_type k = 0; k < d; ++k ) {
				if ( dims2diag_->count( k ) ) {
					(*acc_block_shapes_)({i, k}) = (*acc_block_shapes_)({i - 1, k}) + blocks_->at(i)->shape()[k] ;
				}
				else {
					assert( blocks_->at(i)->shape()[k] == blocks_->at(0)->shape()[k] ) ;
					(*acc_block_shapes_)({i, k}) = blocks_->at(i)->shape()[k] ;
				}
			}
			ndof_ += blocks_->at(i)->size() ;
		}
		shape_ = boost::make_shared<shape_type>( (*acc_block_shapes_)({ndims_type(blocks_->size()) - 1, range<>(0, d)}) ) ;
		dim_factor_ = boost::make_shared<shape_type>( 1, d ) ;
		for ( ndims_type k = 0; k < d; ++k ) {
			size_ *= (*shape_)[k] ;
			if ( k > 0 ) {
				(*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ;
			}
		}
	}

	template <typename = void>
	blocks_array ( std::vector<boost::shared_ptr<array_wrapper<value_type>>> const& blocks, std::set<ndims_type> dims2diag )
		: blocks_( boost::make_shared<std::vector<boost::shared_ptr<array_wrapper<value_type>>>>( blocks ) )
		, dims2diag_( boost::make_shared<std::set<ndims_type>>( std::move( dims2diag ) ) )
		, acc_block_shapes_( )
		, shape_( )
		, dim_factor_( )
		, size_( 1 )
		, ndof_( 0 )
		{
			assert( dims2diag_->size() >= 1 ) ;
			ndims_type d = blocks_->at(0)->shape().size() ;
	    	if ( ndims >= 0 ) {
	    		assert ( ndims == d ) ;
	    	}
			ndof_ += blocks_->at(0)->size() ;
			acc_block_shapes_ = boost::make_shared<dense_matrix<size_type>>( dense_matrix<size_type>( no_init(), {ndims_type(blocks_->size()), d} ) ) ;
			(*acc_block_shapes_)({0, range<>(0, d)}) = blocks_->at(0)->shape() ;
			for ( size_type i = 1; i < blocks_->size(); ++i ) {
				assert( blocks_->at( i )->shape().size() == d ) ;
				for ( ndims_type k = 0; k < d; ++k ) {
					if ( dims2diag_->count( k ) ) {
						(*acc_block_shapes_)({i, k}) = (*acc_block_shapes_)({i - 1, k}) + blocks_->at(i)->shape()[k] ;
					}
					else {
						assert( blocks_->at(i)->shape()[k] == blocks_->at(0)->shape()[k] ) ;
						(*acc_block_shapes_)({i, k}) = blocks_->at(i)->shape()[k] ;
					}
				}
				ndof_ += blocks_->at(i)->size() ;
			}
			shape_ = boost::make_shared<shape_type>( (*acc_block_shapes_)({ndims_type(blocks_->size()) - 1, range<>(0, d)}) ) ;
			dim_factor_ = boost::make_shared<shape_type>( 1, d ) ;
			for ( ndims_type k = 0; k < d; ++k ) {
				size_ *= (*shape_)[k] ;
				if ( k > 0 ) {
					(*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ;
				}
			}
		}

public:
	// Copy constructor -> deep copy
	blocks_array ( blocks_array const& that ) = delete ;

	// Move constructor
	blocks_array ( blocks_array && that ) = default ;

public:
	// Copy assignment -> deep copy
	blocks_array& operator= ( blocks_array const& that ) = delete ;

	// Move assignment
	blocks_array& operator= ( blocks_array&& that ) = delete ;

public:
	// Constructor used in shallow_copy
	blocks_array ( boost::shared_ptr<std::vector<boost::shared_ptr<array_wrapper<value_type>>>> blocks,
			       boost::shared_ptr<std::set<ndims_type>> dims2diag,
			       boost::shared_ptr<dense_matrix<size_type>> acc_block_shapes,
			       boost::shared_ptr<shape_type> shape, boost::shared_ptr<shape_type> dim_factor, size_type size, size_type ndof )
	: blocks_( blocks )
    , dims2diag_( dims2diag )
    , acc_block_shapes_( acc_block_shapes )
    , shape_( shape )
    , dim_factor_( dim_factor )
	, size_( size )
    , ndof_( ndof )
	{}

private:
	boost::shared_ptr<std::vector<boost::shared_ptr<array_wrapper<value_type>>>>    blocks_ ;
	boost::shared_ptr<std::set<ndims_type>>                                         dims2diag_ ;
	boost::shared_ptr<dense_matrix<size_type>>                                      acc_block_shapes_ ;
	boost::shared_ptr<shape_type>                                                   shape_ ;
	boost::shared_ptr<shape_type>                                                   dim_factor_ ;
	size_type                                                                       size_ ;
	size_type                                                                       ndof_ ;

public:
	blocks_array shallow_copy () const {
		return blocks_array( blocks_, dims2diag_, acc_block_shapes_, shape_, dim_factor_, size_, ndof_ );
	}

public:
	shape_type shape () const {
		return *shape_ ;
	}

	size_type size () const {
		return size_ ;
	}

	size_type ndof () const {
		return ndof_ ;
	}

public:
	template <typename I>
	typename std::enable_if< std::is_integral<I>::value, value_type >::type operator[] ( I const& ii ) const {
		assert( ii >= 0 && ii < size_ ) ;
		shape_type j( no_init(), shape_->size() ) ;
		index2multi_index( *dim_factor_, ii, j ) ;
		shape_type j_block( j ) ;
		size_type block_id = -1 ;
		for ( auto k : *dims2diag_ ) {
			if ( block_id == -1 ) {
				block_id = bisect_right( (*acc_block_shapes_)({ range<>(0, blocks_->size()), k }) , j[k] ) ;
			}
			else {
				if ( block_id != bisect_right( (*acc_block_shapes_)({range<>(0, blocks_->size()), k}) , j[k] ) ) {
					return 0 ;
				}
			}
			if ( block_id == 0 ) {
				j_block[k] = j[k] ;
			}
			else {
				j_block[k] = j[k] - (*acc_block_shapes_)( { block_id - 1, k} ) ;
			}
		}
		return (*blocks_->at(block_id))( j_block ) ;
	}

	template <typename X>
	auto operator[] ( X const& x ) const -> typename std::enable_if< is<Array, X>::value, decltype( linear_index_selection( *this, x ) ) >::type {
		return linear_index_selection( *this, x ) ;
	}

	template <typename I>
	auto operator[] ( std::initializer_list<I> const& x ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( linear_index_selection( *this, x ) ) >::type {
		return linear_index_selection( *this, x ) ;
	}

	template <typename J>
	typename std::enable_if< is<Vector, J>::value, value_type >::type operator() ( J const& j ) const {
		assert( shape_->size() == j.size() ) ;
		for ( ndims_type k = 0; k < shape_->size(); ++k ) {
			assert( j[k] >= 0 && j[k] < (*shape_)[k] ) ;
		}
		shape_type j_block( j ) ;
		size_type block_id = -1 ;
		for ( auto k : *dims2diag_ ) {
			if ( block_id == -1 ) {
				block_id = bisect_right( (*acc_block_shapes_)({ range<>(0, blocks_->size()), k }) , j[k] ) ;
			}
			else {
				if ( block_id != bisect_right( (*acc_block_shapes_)({range<>(0, blocks_->size()), k}) , j[k] ) ) {
					return 0 ;
				}
			}
			if ( block_id == 0 ) {
				j_block[k] = j[k] ;
			}
			else {
				j_block[k] = j[k] - (*acc_block_shapes_)( { block_id - 1, k} ) ;
			}
		}
		return (*blocks_->at(block_id))( j_block ) ;
	}

	value_type operator() ( std::initializer_list<size_type> const& j ) const {
		return (*this)( dense_vector<size_type>( j ) ) ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

} ;

template <typename V, std::ptrdiff_t ndims>
struct concept<blocks_array<V, ndims>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseArray
{} ;

template <typename V>
struct concept< blocks_array<V, 0>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseScalar
{} ;

template <typename V>
struct concept< blocks_array<V, 1>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseVector
{} ;

template <typename V>
struct concept< blocks_array<V, 2>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseMatrix
{} ;

} // namespace glas3


#endif
