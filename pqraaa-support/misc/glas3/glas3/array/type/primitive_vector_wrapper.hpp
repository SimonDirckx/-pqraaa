//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_primitive_vector_wrapper_hpp
#define glas3_array_type_primitive_vector_wrapper_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/primitive_scalar_wrapper.hpp>

#include <functional>
#include <initializer_list>
#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <utility>
#include <iostream>

namespace glas3 {

template <typename V>
class primitive_vector_wrapper {
	static_assert( std::is_arithmetic<V>::value || boost::is_complex<V>::value, "V should be arithmetic" ) ;

public:
	typedef V                                           value_type ;
	typedef primitive_scalar_wrapper<std::ptrdiff_t>    shape_type ;
	typedef typename shape_type::value_type             size_type ;
	typedef typename shape_type::size_type              ndims_type ;

public:
	explicit primitive_vector_wrapper ()
	: array_( boost::make_shared< boost::any >( ) )
	, shape_( boost::make_shared< shape_type >( 0 ) )
	, size_( 0 )
	, ndof_( 0 )
	{}

	primitive_vector_wrapper ( value_type const& value )
	: array_( )
	, shape_( boost::make_shared< shape_type >( 1 ) )
	, size_( 1 )
	, ndof_( 1 )
	{
		auto array = boost::make_shared<value_type[]>( 1 ) ;
		array[0] = value ;
		array_ = boost::make_shared< boost::any >( std::move( array ) ) ;
		indexer_ = boost::make_shared<std::function<value_type ( size_type const& )>> (
				[this] ( size_type const& i ) { return ( boost::any_cast<boost::shared_ptr<value_type[]>>( *( this->array_ ) ) ) [0] ; } ) ;
	}

	// Copy constructor from std::initializer_list< value_type > -> deep copy
	primitive_vector_wrapper ( std::initializer_list< value_type > const& data )
	: array_( )
	, shape_( boost::make_shared< shape_type >( data.size() ) )
	, size_( data.size() )
	, ndof_( data.size() )
	{
		auto array = boost::make_shared<value_type[]>( data.size() ) ;
		size_type i = 0 ;
		for ( auto x = data.begin(); x != data.end(); ++x, ++i) {
			array[i] = *x ;
		}
		array_ = boost::make_shared< boost::any >( std::move( array ) ) ;
		indexer_ = boost::make_shared<std::function<value_type ( size_type const& )>> (
				[this] ( size_type const& i ) { return ( boost::any_cast<boost::shared_ptr<value_type[]>>( *( this->array_ ) ) )[ i ] ; } ) ;
	}

	// Copy constructor from Array -> shallow copy
	template <typename E>
	primitive_vector_wrapper ( E const& that, typename std::enable_if<is<Array, E>::value>::type* = 0 )
	: array_( boost::make_shared< boost::any >( boost::make_shared<E>( std::move( that.shallow_copy() ) ) ) )
	  , shape_( boost::make_shared< shape_type >( that.size() ) )
	  , size_( that.size() )
	  , ndof_( that.ndof() )
	  {
		indexer_ = boost::make_shared<std::function<value_type ( size_type const& )>> (
				[this] ( size_type const& i ) { return ( *boost::any_cast<boost::shared_ptr<E>>( *( this->array_ ) ) ) [ i ]  ; } ) ;
	  }

public:
	// Copy constructor -> shallow copy
	primitive_vector_wrapper ( primitive_vector_wrapper const& that ) = default ;

	// Move constructor
	primitive_vector_wrapper ( primitive_vector_wrapper && ) = default ;

public:
	// Copy assignment -> do not allow
	primitive_vector_wrapper& operator= ( primitive_vector_wrapper const& that ) = delete ; // not assignable

	// Move assignment
	primitive_vector_wrapper& operator= ( primitive_vector_wrapper&& that ) = delete ;

public:
	// Constructor used in shallow_copy
	primitive_vector_wrapper ( boost::shared_ptr<std::function<value_type ( size_type const& )>> indexer, boost::shared_ptr<boost::any> array, boost::shared_ptr<shape_type> shape, size_type size, size_type ndof )
	: indexer_( indexer )
	, array_( array )
	, shape_( shape )
	, size_( size )
	, ndof_( ndof )
	{}

private:
	boost::shared_ptr<std::function<value_type ( size_type const& )>>          indexer_;
	boost::shared_ptr<boost::any>                                              array_ ;
	boost::shared_ptr<shape_type>                                              shape_ ;
	size_type                                                                  size_ ;
	size_type                                                                  ndof_ ;

public:
	primitive_vector_wrapper shallow_copy () const {
		return primitive_vector_wrapper( indexer_, array_, shape_, size_, ndof_ );
	}

public:
	shape_type const& shape () const {
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
	typename std::enable_if< std::is_integral<I>::value, value_type>::type operator[] ( I const& i ) const {
		assert( i >= 0 && i < size_ ) ;
		return (*indexer_)( i ) ;
	}

	template <typename X>
	typename std::enable_if< is<Vector, X>::value, value_type >::type operator() ( X const& x ) const {
		assert( x.size() == 1 ) ;
		assert( x[0] >= 0 && x[0] < size_ ) ;
		return (*indexer_)( x[0] ) ;
	}

	value_type operator() ( std::initializer_list<size_type> const& x ) const {
		assert( x.size() == 1 ) ;
		assert( *x.begin() >= 0 && *x.begin() < size_ ) ;
		return (*indexer_)( *x.begin() ) ;
	}

} ;

template <typename V>
struct concept<primitive_vector_wrapper<V>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseVector
  {} ;

} ; // namespace glas3

#endif
