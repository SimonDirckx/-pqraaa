//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_indexed_array_hpp
#define glas3_array_type_indexed_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/dense_array/concept/index.hpp>

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

namespace glas3 {

template <typename T, typename S>
class indexed_array {
	static_assert( is<Array, T>::value, "T should be an Array" ) ;
	static_assert( is<Index, S>::value, "S should be an Index" ) ;

public:
	typedef T                               array_type ;
	typedef S                               index_type ;
	typedef typename T::value_type          value_type ;
	typedef typename S::shape_type          shape_type ;
	typedef typename S::size_type           size_type ;
	typedef typename S::ndims_type          ndims_type ;

public:
	indexed_array ( )
	: array_( boost::make_shared< array_type >( ) )
	, index_( boost::make_shared< index_type >( ) )
	{}

	indexed_array ( array_type const& array, index_type const& index )
	: array_( boost::make_shared< array_type >( std::move( array.shallow_copy() ) ) )
	, index_( boost::make_shared< index_type >( std::move( index.shallow_copy() ) ) )
	{}

public:
	// Copy constructor -> deep copy
	indexed_array ( indexed_array const& that )
    : array_( boost::make_shared< array_type >( *that.array_ ) )
	, index_( boost::make_shared< index_type >( *that.index_ ) )
    {}

	// Move constructor
	indexed_array ( indexed_array && that ) = default;

public:
	// Copy assignment -> deep copy
	indexed_array& operator= ( indexed_array const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	// Move assignment
	indexed_array& operator= ( indexed_array&& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	indexed_array& operator= ( value_type const& value ) {
		assert( size() == 1 );
		(*this)[0] = value ;
		return *this;
	}

	indexed_array& operator= ( std::initializer_list<value_type> const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, indexed_array& >::type operator= ( E const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, indexed_array& >::type operator= ( E const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

public:
	// Constructor used in shallow_copy
	indexed_array ( boost::shared_ptr<array_type> array, boost::shared_ptr<index_type> index )
: array_( array )
, index_( index )
{}

private:
	boost::shared_ptr<array_type>     array_ ;
	boost::shared_ptr<index_type>     index_ ;

public:
	indexed_array shallow_copy () const {
		return indexed_array( array_, index_ );
	}

public:
	auto shape () const -> decltype( index_->shape_in() ) {
		return index_->shape_in() ;
	}

	auto size () const -> decltype( index_->lin_size() ) {
		return index_->lin_size() ;
	}

	auto ndof () const -> decltype( array_->ndof() + index_->ndof() ) {
		return array_->ndof() + index_->ndof() ;
	}

public:
	template <typename I>
	auto operator[] ( I const& i ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( (*array_)(*index_) ) >::type {
		index_->set_lin_in( i ) ;
		return (*array_)(*index_) ;
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
	auto operator() ( J const& j ) const -> typename std::enable_if< is<Vector, J>::value, decltype ( (*array_)(*index_) ) >::type {
		index_->set_multi_in( j ) ;
		return (*array_)(*index_) ;
	}

	auto operator() ( std::initializer_list<size_type> const& j ) const -> decltype ( (*array_)(*index_) ) {
		index_->set_multi_in( j ) ;
		return (*array_)(*index_) ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

} ;

template <typename T, typename S>
struct concept<indexed_array<T, S>, typename std::enable_if< is<Array, T>::value && is<Index, S>::value >::type >
: DenseArray
  {};

} // namespace glas3

#endif
