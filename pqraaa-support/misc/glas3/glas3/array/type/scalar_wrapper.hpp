//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_scalar_wrapper_hpp
#define glas3_array_type_scalar_wrapper_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/empty_array.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>

#include <functional>
#include <initializer_list>
#include <type_traits>
#include <cassert>
#include <vector>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <utility>

namespace glas3 {

template <typename V>
class scalar_wrapper {
	static_assert( std::is_arithmetic<V>::value || boost::is_complex<V>::value, "V should be arithmetic" ) ;

public:
	typedef V                                 value_type ;
	typedef empty_array                       shape_type ;
	typedef typename shape_type::value_type   size_type ;
	typedef typename shape_type::size_type    ndims_type ;

public:
	explicit scalar_wrapper ()
	: array_( boost::make_shared< boost::any >( ) )
	, shape_( boost::make_shared< shape_type >( ) )
	, size_( 0 )
	, ndof_( 0 )
	{}

	scalar_wrapper ( value_type const& value )
	: array_( boost::make_shared< boost::any >( value ) )
	, shape_( boost::make_shared< shape_type >( ) )
	, size_( 1 )
	, ndof_( 1 )
	{
		indexer_ = boost::make_shared<std::function<value_type ( size_type const& )>> (
				[this] ( size_type const& i ) { return ( boost::any_cast<value_type>( *( this->array_ ))) ; } ) ;
	}

	// Copy constructor from std::initializer_list< value_type > -> deep copy
	scalar_wrapper ( std::initializer_list< value_type > const& data )
	: array_( )
	, shape_( boost::make_shared< shape_type >( ) )
	, size_( data.size() )
	, ndof_( data.size() )
	{
		assert( data.size() <= 1 ) ;
		if ( size_ == 1 ) {
			array_ = boost::make_shared< boost::any >( std::move( value_type( *data.begin() ) ) ) ;
			indexer_ = boost::make_shared<std::function<value_type ( size_type const& )>> (
					[this] ( size_type const& i ) { return ( boost::any_cast<value_type>( *( this->array_ ))) ; } ) ;
		}
	}

	//        // Copy constructor from Container -> deep copy
	//        template <typename E>
	//        scalar_wrapper ( E const& data, typename std::enable_if<is<Container, E>::value>::type* = 0 )
	//        : array_( )
	//        , shape_( boost::make_shared< shape_type >( ) )
	//        , size_( data.size() )
	//        , ndof_( data.size() )
	//        {
	//        	assert( data.size() <= 1 ) ;
	//        	if ( size_ == 1 ) {
	//        	    array_ = boost::make_shared< boost::any >( std::move( value_type( *data.begin() ) ) ) ;
	//        	    indexer_ = [this] ( size_type const& i ) { return ( boost::any_cast<value_type>( *( this->array_ ))) ; } ;
	//            }
	//        }
	//
	//        // Move constructor from RandomAccessContainer -> move
	//        template <typename E>
	//        scalar_wrapper ( E && data, typename std::enable_if<is<RandomAccessContainer, E>::value>::type* = 0 )
	//        : array_( )
	//        , shape_( boost::make_shared< shape_type >( ) )
	//        , size_( data.size() )
	//        , ndof_( data.size() )
	//        {
	//        	assert( data.size() <= 1 ) ;
	//            if ( size_ == 1 ) {
	//        	    array_ = boost::make_shared< boost::any >( std::move( data ) ) ;
	//        	    indexer_ = [this] ( size_type const& i ) { return ( boost::any_cast<E>( *( this->array_ ))) [ i ]  ; } ;
	//            }
	//        }

	// Copy constructor from Array -> shallow copy
	template <typename E>
	scalar_wrapper ( E const& that, typename std::enable_if<is<Array, E>::value>::type* = 0 )
	: array_( boost::make_shared< boost::any >( std::move( that.shallow_copy() ) ) )
	  , shape_( boost::make_shared< shape_type >( ) )
	  , size_( that.size() )
	  , ndof_( that.ndof() )
	  {
		assert( that.size() <= 1 ) ;
		indexer_ = boost::make_shared<std::function<value_type ( size_type const& )>> (
				[this] ( size_type const& i ) { return ( boost::any_cast<E>( *( this->array_ ))) [ i ]  ; } ) ;
	  }

public:
	// Copy constructor -> shallow copy
	scalar_wrapper ( scalar_wrapper const& that ) = default ;

	// Move constructor
	scalar_wrapper ( scalar_wrapper && ) = default ;

public:
	// Copy assignment -> do not allow
	scalar_wrapper& operator= ( scalar_wrapper const& that ) = delete ; // not assignable

	// Move assignment
	scalar_wrapper& operator= ( scalar_wrapper&& that ) = delete ;

public:
	// Constructor used in shallow_copy
	scalar_wrapper ( boost::shared_ptr<std::function<value_type ( size_type const& )>> indexer, boost::shared_ptr<boost::any> array, boost::shared_ptr<shape_type> shape, size_type size, size_type ndof )
	: indexer_( indexer )
	, array_( array )
	, shape_( shape )
	, size_( size )
	, ndof_( ndof )
	{}

private:
	boost::shared_ptr<std::function<value_type ( size_type const& )>>     indexer_;
	boost::shared_ptr<boost::any>                                         array_ ;
	boost::shared_ptr<shape_type>                                         shape_ ;
	size_type                                                             size_ ;
	size_type                                                             ndof_ ;

public:
	scalar_wrapper shallow_copy () const {
		return scalar_wrapper( indexer_, array_, shape_, size_, ndof_ );
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
	auto operator[] ( X const& x ) const -> typename std::enable_if< is<Array, X>::value, decltype( linear_index_selection( *this, x ) ) >::type {
		return linear_index_selection( *this, x ) ;
	}

	template <typename I>
	auto operator[] ( std::initializer_list<I> const& x ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( linear_index_selection( *this, x ) ) >::type {
		return linear_index_selection( *this, x ) ;
	}

	template <typename J>
	typename std::enable_if< is<Vector, J>::value, value_type >::type operator() ( J const& j ) const {
		assert( j.size() == 0 ) ;
		assert( 0 < size_ ) ;
		return (*indexer_)( 0 ) ;
	}

	value_type operator() ( std::initializer_list<size_type> const& j ) const {
		assert( j.size() == 0 ) ;
		assert( 0 < size_ ) ;
		return (*indexer_)( 0 ) ;
	}

} ;

template <typename V>
struct concept<scalar_wrapper<V>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseScalar
  {};

} ; // namespace glas3

#endif
