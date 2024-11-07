//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_glas2_dense_matrix_hpp
#define glas3_array_dense_array_type_glas2_dense_matrix_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas2/concept/is.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>

#include <glas3/array/type/primitive_vector_wrapper.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection.hpp>
#include <glas3/array/dense_array/algorithm/assign.hpp>

#include <vector>
#include <cassert>
#include <type_traits>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>
#include <iostream>
#include <initializer_list>

namespace glas3 {

template <typename T>
class glas2_dense_matrix {
public:
	typedef T                                        array_type ;
	typedef typename T::value_type                   value_type ;
	typedef dense_vector<typename T::size_type>      shape_type ;
	typedef typename shape_type::value_type          size_type ;
	typedef typename shape_type::size_type           ndims_type ;

public:
	glas2_dense_matrix ()
    : array_( )
    , shape_( boost::make_shared< shape_type >( std::move ( shape_type ( { 0, 0 } ) ) ) )
    , size_( 0 )
    {}

	glas2_dense_matrix ( array_type& array )
	: array_( array )
    , shape_( boost::make_shared< shape_type >( std::move ( shape_type ( { size_type( array.num_rows() ), size_type( array.num_columns() ) } ) ) ) )
    , size_( array.num_rows() * array.num_columns() )
	{}

public:
	// Copy constructor -> do not allow
	glas2_dense_matrix ( glas2_dense_matrix const& that ) = delete ;

	// Move constructor
	glas2_dense_matrix ( glas2_dense_matrix && that ) = default ;

public:
	// Copy assignment -> deep copy
	glas2_dense_matrix& operator= ( glas2_dense_matrix const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	// Move assignment
	glas2_dense_matrix& operator= ( glas2_dense_matrix&& that ) = default;

	glas2_dense_matrix& operator= ( value_type const& value ) {
		assert( size() == 1 );
		array_( 0, 0 ) = value ;
		return *this;
	}

	glas2_dense_matrix& operator= ( std::initializer_list<value_type> const& that ) {
		assert( that.size() == size() );
		auto that_it = that.begin() ;
		for ( size_type j = 0; j < array_.num_columns(); ++j ) {
			for ( size_type i = 0; i < array_.num_rows(); ++i, ++that_it ) {
				array_( i, j ) = *that_it ;
			}
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, glas2_dense_matrix& >::type operator= ( E const& that ) {
		assert( that.size() == size() );
		auto that_it = that.begin() ;
		for ( size_type j = 0; j < array_.num_columns(); ++j ) {
			for ( size_type i = 0; i < array_.num_rows(); ++i, ++that_it ) {
				array_( i, j ) = *that_it ;
			}
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, glas2_dense_matrix& >::type operator= ( E const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

public:
	// Constructor used in shallow_copy
	glas2_dense_matrix ( array_type& array, boost::shared_ptr<shape_type> shape, size_type size )
    : array_( array )
    , shape_( shape )
    , size_( size )
    {}

private:
	array_type&                      array_ ;
	boost::shared_ptr<shape_type>    shape_ ;
	size_type                        size_ ;

public:
	glas2_dense_matrix shallow_copy () const {
		return glas2_dense_matrix( array_, shape_, size_ );
	}

public:
	shape_type const& shape () const {
		return *shape_;
	}

	size_type size () const {
		return size_;
	}

	size_type ndof () const {
		return size_;
	}

public:
    template <typename I>
	auto operator[] ( I const& i ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( array_( size_type(), size_type() ) ) >::type {
		assert( i >= 0 && i < size_ ) ;
	   	std::ldiv_t dv = std::ldiv( i, array_.num_rows() ) ;
		return array_( dv.rem, dv.quot ) ;
	}

    template <typename X>
    auto operator[] ( X const& x ) const -> typename std::enable_if< is<Array, X>::value, decltype ( linear_index_selection( *this, x ) ) >::type {
    	return linear_index_selection( *this, x ) ;
    }

    template <typename I>
    auto operator[] ( std::initializer_list<I> const& x ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( linear_index_selection( *this, x ) ) >::type {
      	return linear_index_selection( *this, x ) ;
    }

    template <typename X>
    auto operator() ( X const& x ) const -> typename std::enable_if< is<Vector, X>::value, decltype( array_( size_type(), size_type() ) ) >::type {
    	assert( x.size() == 2 ) ;
    	return array_( x[0], x[1] ) ;
    }

	auto operator() ( std::initializer_list<size_type> const& x ) const -> decltype( array_( size_type(), size_type() ) ) {
		assert( x.size() == 2 ) ;
		return array_( *x.begin(), *(++x.begin()) ) ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

} ;

template <typename T>
struct concept< glas2_dense_matrix<T>, typename std::enable_if< glas2::is< glas2::DenseMatrix, T >::value >::type >
: DenseMatrix
{};

} // namespace glas3


#endif
