//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_indirect_array_hpp
#define glas3_array_type_indirect_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/primitive_vector_wrapper.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection.hpp>
#include <glas3/array/dense_array/algorithm/assign.hpp>

#include <initializer_list>
#include <vector>
#include <cassert>
#include <type_traits>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>
#include <iostream>

namespace glas3 {

template <typename T, typename S>
class indirect_array {
	static_assert( is<Array, T>::value, "T should be an Array" ) ;
	static_assert( is<Array, S>::value, "S should be an Array" ) ;

public:
	typedef T                               array_type ;
	typedef S                               indexer_type ;
	typedef typename T::value_type          value_type ;
	typedef typename S::shape_type          shape_type ;
	typedef typename S::size_type           size_type ;
	typedef typename S::ndims_type          ndims_type ;

public:
	indirect_array ( )
	: array_( boost::make_shared< array_type >( ) )
	, indexer_( boost::make_shared< indexer_type >( ) )
	{}

	indirect_array ( array_type const& array, indexer_type const& indexer )
	: array_( boost::make_shared< array_type >( std::move( array.shallow_copy() ) ) )
	, indexer_( boost::make_shared< indexer_type >( std::move( indexer.shallow_copy() ) ) )
	{}

public:
	// Copy constructor -> deep copy
	indirect_array ( indirect_array const& that )
    : array_( boost::make_shared< array_type >( *that.array_ ) )
    , indexer_( boost::make_shared< indexer_type >( *that.indexer_ ) )
    {}

	// Move constructor
	indirect_array ( indirect_array && that ) = default;

public:
	// Copy assignment -> deep copy
	indirect_array& operator= ( indirect_array const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	// Move assignment
	indirect_array& operator= ( indirect_array&& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	indirect_array& operator= ( value_type const& value ) {
		assert( size() == 1 );
		(*this)[0] = value ;
		return *this;
	}

	indirect_array& operator= ( std::initializer_list<value_type> const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, indirect_array& >::type operator= ( E const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, indirect_array& >::type operator= ( E const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

public:
	// Constructor used in shallow_copy
	indirect_array ( boost::shared_ptr<array_type> array, boost::shared_ptr<indexer_type> indexer )
	: array_( array )
	, indexer_( indexer )
	{}

private:
	boost::shared_ptr<array_type>     array_ ;
	boost::shared_ptr<indexer_type>   indexer_ ;

public:
	indirect_array shallow_copy () const {
		return indirect_array( array_, indexer_ );
	}

public:
	auto shape () const -> decltype( indexer_->shape() ) {
		return indexer_->shape() ;
	}

	auto size () const -> decltype( indexer_->size() ) {
		return indexer_->size() ;
	}

	auto ndof () const -> decltype( array_->ndof() + indexer_->ndof() ) {
		return array_->ndof() + indexer_->ndof() ;
	}

public:
	template <typename I>
	auto operator[] ( I const& i ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( (*array_)[(*indexer_)[i]] ) >::type {
		return (*array_)[(*indexer_)[i]] ;
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
	auto operator() ( J const& j ) const -> typename std::enable_if< is<Vector, J>::value, decltype ( (*array_)[(*indexer_)( j )] ) >::type {
		return (*array_)[(*indexer_)( j )] ;
	}

	auto operator() ( std::initializer_list<size_type> const& j ) const -> decltype ( (*array_)[(*indexer_)( j )] ) {
		return (*array_)[(*indexer_)( j )] ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

} ;

template <typename T, typename S>
struct concept<indirect_array<T, S>, typename std::enable_if< is<DenseArray, T>::value && is<DenseArray, S>::value && !is<DenseMatrix, S>::value && !is<DenseVector, S>::value && !is<DenseScalar, S>::value >::type >
: DenseArray
  {};

template <typename T, typename S>
struct concept<indirect_array<T, S>, typename std::enable_if< is<DenseArray, T>::value && is<DenseMatrix, S>::value >::type >
: DenseMatrix
  {};

template <typename T, typename S>
struct concept<indirect_array<T, S>, typename std::enable_if< is<DenseArray, T>::value && is<DenseVector, S>::value >::type >
: DenseVector
  {};

template <typename T, typename S>
struct concept<indirect_array<T, S>, typename std::enable_if< is<DenseArray, T>::value && is<DenseScalar, S>::value >::type >
: DenseScalar
  {};

} // namespace glas3


#endif
