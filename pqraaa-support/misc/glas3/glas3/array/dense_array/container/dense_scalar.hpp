//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_container_dense_scalar_hpp
#define glas3_array_dense_array_container_dense_scalar_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/contiguous_dense_array.hpp>

#include <glas3/array/type/empty_array.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>

#include <initializer_list>
#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <vector>
#include <utility>

namespace glas3 {

template <typename V>
class dense_scalar {
	static_assert( std::is_arithmetic<V>::value || boost::is_complex<V>::value, "V should be arithmetic" ) ;

public:
	typedef V                value_type ;
	typedef empty_array      shape_type ;
	typedef std::ptrdiff_t   size_type ;
	typedef std::ptrdiff_t   ndims_type ;

public:
	dense_scalar ()
    : data_( boost::make_shared_noinit< value_type[] > ( 0 ) )
    , shape_( boost::make_shared< empty_array > ( ) )
    , size_( 0 )
    {}

	// Copy constructor from value_type -> deep copy
	dense_scalar ( value_type const& value )
	: data_( boost::make_shared< value_type[] > ( 1, value ) )
	, shape_( boost::make_shared< empty_array > ( ) )
	, size_( 1 )
	{}

	// Copy constructor from std::initializer_list< value_type > -> deep copy
	dense_scalar ( std::initializer_list< value_type > const& that )
	: data_( boost::make_shared_noinit< value_type[] >( that.size() ) )
	, shape_( boost::make_shared< empty_array >( ) )
	, size_( that.size() )
	{
		assert( that.size() <= 1 ) ;
		if ( size_ == 1 ) { data_[0] = *that.begin() ; }
	}

	// Copy constructor from Container -> deep copy
	template <typename E>
	dense_scalar ( E const& that, typename std::enable_if<is<Container, E>::value>::type* = 0 )
	: data_( boost::make_shared_noinit< value_type[] >( that.size() ) )
	  , shape_( boost::make_shared< empty_array >( ) )
	  , size_( that.size() )
	  {
		assert( that.size() <= 1 ) ;
		if ( size_ == 1 ) { data_[0] = *that.begin() ; }
	  }

	// Move constructor from std::vector< value_type > -> move ?? how ??
	//        dense_scalar ( std::vector< value_type > && that )
	//        : data_( boost::make_shared_noinit< value_type[] >( that.size() ) )
	//        , shape_( boost::make_shared< empty_array >( ) )
	//        , size_( that.size() )
	//        {
	//          	assert( that.size() <= 1 ) ;
	//          	if ( size_ == 1 ) { data_[0] = that[0] ; }
	//        }

	// Copy constructor from Array -> deep copy
	template <typename E>
	dense_scalar ( E const& that, typename std::enable_if<is<Array, E>::value>::type* = 0 )
	: data_( boost::make_shared_noinit< value_type[] > ( that.size() ) )
	  , shape_( boost::make_shared< empty_array > ( ) )
	  , size_( that.size() )
	  {
		assert( that.size() <= 1 ) ;
		if ( size_ == 1 ) { data_[0] = that[0] ; }
	  }

public:
	// Copy constructor -> deep copy
	dense_scalar ( dense_scalar const& that )
    : data_( boost::make_shared_noinit< value_type[] > ( that.size_ ) )
    , shape_( boost::make_shared< empty_array > ( ) )
    , size_( that.size_ )
    {
		if ( size_ == 1 ) { data_[0] = that.data_[0] ; }
    }

	// Move constructor
	dense_scalar ( dense_scalar && ) = default ;

public:
	// Copy assignment -> deep copy
	dense_scalar& operator= ( dense_scalar const& that ) {
		assert( that.size_ == size_ ) ;
		if ( size_ == 1 ) { data_[0] = that.data_[0] ; }
		return *this ;
	}

	// Move assignment
	dense_scalar& operator= ( dense_scalar&& that ) {
		assert( that.size_ == size_ ) ;
		if ( size_ == 1 ) { data_[0] = that.data_[0] ; }
		return *this ;
	}

	dense_scalar& operator= ( value_type const& value ) {
		assert( size_ == 1 ) ;
		data_[0] = value ;
		return *this;
	}

	dense_scalar& operator= ( std::initializer_list< value_type > const& that ) {
		assert( that.size() == size_ );
		if ( size_ == 1 ) { data_[0] = *that.begin() ; }
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, dense_scalar& >::type operator= ( E const& that ) {
		assert( that.size() == size_ );
		if ( size_ == 1 ) { data_[0] = *that.begin() ; }
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, dense_scalar& >::type operator= ( E const& that ) {
		assert( that.size() == size_ ) ;
		if ( size_ == 1 ) { data_[0] = that[0] ; }
		return *this ;
	}

private:
	// Constructor used in shallow_copy
	dense_scalar ( boost::shared_ptr< value_type[] > data, boost::shared_ptr< empty_array > shape, size_type size )
    : data_( data )
    , shape_( shape )
    , size_( size )
    {}

private:
	boost::shared_ptr< value_type[] >  data_ ;
	boost::shared_ptr< empty_array >   shape_ ;
	size_type                          size_ ;

public:
	dense_scalar shallow_copy () const {
		return dense_scalar( data_, shape_, size_ );
	}

public:
	boost::shared_ptr< value_type[] > const data_ptr () const {
		return data_ ;
	}

public:
	shape_type const& shape () const {
		return *shape_ ;
	}

	size_type size () const {
		return size_ ;
	}

	size_type ndof () const {
		return size_ ;
	}

public:
	template <typename I>
	typename std::enable_if< std::is_integral<I>::value, value_type&>::type operator[] ( I const& i ) const {
		assert( i >= 0 && i < size_ ) ;
		return data_[0] ;
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
	typename std::enable_if< is<Vector, J>::value, value_type& >::type operator() ( J const& j ) const {
		assert( j.size() == 0 ) ;
		assert( 0 < size_ ) ;
		return data_[0] ;
	}

	value_type& operator() ( std::initializer_list<size_type> const& j ) const {
		assert( j.size() == 0 ) ;
		assert( 0 < size_ ) ;
		return data_[0] ;
	}

	//    	template < typename S >
	//    	dense_scalar operator() ( std::initializer_list<S> const& j ) const {
		//    		assert( j.size() == 0 ) ;
		//        	return this->shallow_copy() ;
		//    	}

} ;

template <typename V>
struct concept< dense_scalar<V>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: ContiguousDenseScalar
  {};

} ; // namespace glas3

#endif
