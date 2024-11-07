//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_glas2_dense_vector_hpp
#define glas3_array_dense_array_type_glas2_dense_vector_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>

#include <glas3/array/type/primitive_vector_wrapper.hpp>
#include <glas3/array/dense_array/container/dense_scalar.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection_of_vector.hpp>
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

template <typename T>
class glas2_dense_vector {
public:
	typedef T                                        array_type ;
	typedef typename T::value_type                   value_type ;
	typedef dense_scalar<typename T::size_type>      shape_type ;
	typedef typename shape_type::value_type          size_type ;
	typedef typename shape_type::size_type           ndims_type ;

public:
	glas2_dense_vector ()
    : array_( )
    , shape_( boost::make_shared< shape_type >( 0 ) )
    , size_( 0 )
    {}

	glas2_dense_vector ( array_type& array )
	: array_( array )
    , shape_( boost::make_shared< shape_type >( array.size() ) )
    , size_( array.size() )
	{}

public:
	// Copy constructor -> do not allow
	glas2_dense_vector ( glas2_dense_vector const& that ) = delete ;

	// Move constructor
	glas2_dense_vector ( glas2_dense_vector && that ) = default ;

public:
	// Copy assignment -> deep copy
	glas2_dense_vector& operator= ( glas2_dense_vector const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	// Move assignment
	glas2_dense_vector& operator= ( glas2_dense_vector&& that ) = default;

	glas2_dense_vector& operator= ( value_type const& value ) {
		assert( size() == 1 );
		array_( 0 ) = value ;
		return *this;
	}

	glas2_dense_vector& operator= ( std::initializer_list<value_type> const& that ) {
		assert( that.size() == size() );
		auto that_it = that.begin() ;
		for ( size_type i = 0; i < array_.size(); ++i, ++that_it ) {
			array_( i ) = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, glas2_dense_vector& >::type operator= ( E const& that ) {
		assert( that.size() == size() );
		auto that_it = that.begin() ;
		for ( size_type i = 0; i < array_.num_rows(); ++i, ++that_it ) {
			array_( i ) = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, glas2_dense_vector& >::type operator= ( E const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

public:
	// Constructor used in shallow_copy
	glas2_dense_vector ( array_type& array, boost::shared_ptr<shape_type> shape, size_type size )
    : array_( array )
    , shape_( shape )
    , size_( size )
    {}

private:
	array_type&                      array_ ;
	boost::shared_ptr<shape_type>    shape_ ;
	size_type                        size_ ;

public:
	glas2_dense_vector shallow_copy () const {
		return glas2_dense_vector( array_, shape_, size_ );
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
	auto operator[] ( I const& i ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( array_( i ) ) >::type {
		assert( i >= 0 && i < size_ ) ;
		return array_( i ) ;
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
    auto operator() ( X const& x ) const -> typename std::enable_if< is<Vector, X>::value, decltype( array_( x[0] ) ) >::type {
    	assert( x.size() == 1 ) ;
    	return array_( x[0] ) ;
    }

	auto operator() ( std::initializer_list<size_type> const& x ) const -> decltype( array_( *x.begin() ) ) {
		assert( x.size() == 1 ) ;
		return array_( *x.begin() ) ;
	}

	template <typename = void>
	auto operator() ( std::vector<primitive_vector_wrapper<size_type>> const& s ) const -> decltype ( block_selection_of_vector( *this, *s.begin() ) ) {
		assert( s.size() == 1 ) ;
		return block_selection_of_vector( *this, *s.begin() ) ;
	}
} ;

template <typename T>
struct concept< glas2_dense_vector<T>, typename std::enable_if< glas2::is< glas2::DenseVector, T >::value >::type >
: DenseVector
{};

} // namespace glas3


#endif
