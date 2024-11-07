//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_container_dense_vector_hpp
#define glas3_array_dense_array_container_dense_vector_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/contiguous_dense_array.hpp>

#include <glas3/array/type/no_init.hpp>
#include <glas3/array/type/primitive_vector_wrapper.hpp>
#include <glas3/array/dense_array/container/dense_scalar.hpp>

#include <glas3/array/dense_array/algorithm/assign.hpp>
#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection_of_vector.hpp>

#include <initializer_list>
#include <vector>
#include <type_traits>
#include <cassert>
#include <algorithm>
#include <utility>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_complex.hpp>

namespace glas3 {

template <typename V>
class dense_vector {
	static_assert( std::is_arithmetic<V>::value || boost::is_complex<V>::value, "V should be arithmetic" ) ;

public:
	typedef V                                 value_type ;
	typedef dense_scalar<std::ptrdiff_t>      shape_type ;
	typedef typename shape_type::value_type   size_type ;
	typedef typename shape_type::size_type    ndims_type ;

public:
	dense_vector ()
    : data_( boost::make_shared_noinit< value_type[] >( 0 ) )
    , shape_( boost::make_shared< shape_type >( 0 ) )
    , size_( 0 )
    {}

	// Copy constructor from value_type -> deep copy
	dense_vector ( value_type const& value )
	: data_( boost::make_shared< value_type[] >( 1, value ) )
	, shape_( boost::make_shared< shape_type >( 1 ) )
	, size_( 1 )
	{}

	// Copy constructor from std::initializer_list< value_type > -> deep copy
	dense_vector ( std::initializer_list< value_type > const& data )
	: data_( boost::make_shared_noinit< value_type[] >( data.size() ) ) // boost::make_shared< value_type[] >( data ) doesn't work ??
	, shape_( boost::make_shared< shape_type >( data.size() ) )
	, size_( data.size() )
	{
		//std::copy( data.begin(), data.end(), data_ ) ; -> std::copy doesn't seem to work with boost::shared_ptr ??
		size_type i = 0 ;
		for ( auto x = data.begin(); x != data.end(); ++x, ++i) {
			data_[i] = *x ;
		}
	}

	// Copy constructor from Container -> deep copy
	template <typename E>
	dense_vector ( E const& that, typename std::enable_if<is<Container, E>::value>::type* = 0 )
	: data_( boost::make_shared_noinit< value_type[] >( that.size() ) )
	  , shape_( boost::make_shared< shape_type >( that.size() ) )
	  , size_( that.size() )
	  {
		size_type i = 0 ;
		for ( auto x = that.begin(); x != that.end(); ++x, ++i) {
			data_[i] = *x ;
		}
	  }

	// uninitialized dense_vector constructor
	dense_vector ( no_init, shape_type const& size )
	: data_( boost::make_shared_noinit< value_type[] >( size[0] ) )
	, shape_( boost::make_shared< shape_type >( size[0] ) )
	, size_( size[0] )
	{}

	dense_vector ( value_type const& value, shape_type const& size )
	: data_( boost::make_shared< value_type[] >( size[0], value ) )
	, shape_( boost::make_shared< shape_type >( size[0] ) )
	, size_( size[0] )
	{}

	// Copy constructor from Array -> deep copy
	template <typename E>
	dense_vector ( E const& that, typename std::enable_if<is<Array, E>::value>::type* = 0 )
	: data_( boost::make_shared_noinit< value_type[] >( that.size() ) )
	  , shape_( boost::make_shared< shape_type >( that.size() ) )
	  , size_( that.size() )
	  {
		assign( *this, that ) ;
	  }

public:
	// Copy constructor -> deep copy
	dense_vector ( dense_vector const& that )
	: data_( boost::make_shared_noinit< value_type[] >( that.size_ ) )
	, shape_( boost::make_shared< shape_type >( that.size_ ) )
	, size_( that.size_ )
	{
			//std::copy( that.data_, that.data_ + that.size_, data_ ) ; -> std::copy doesn't seem to work with boost::shared_ptr ??
			assign( *this, that ) ;
	}

	// Move constructor
	dense_vector ( dense_vector && ) = default ;

public:
	// Copy assignment -> deep copy
	dense_vector& operator= ( dense_vector const& that ) {
		assert( that.size_ == size_ ) ;
		//std::copy( that.data_, that.data_ + that.size_, data_ ) ; -> std::copy doesn't seem to work with boost::shared_ptr ??
		assign( *this, that ) ;
		return *this ;
	}

	// Move assignment
	dense_vector& operator= ( dense_vector&& that ) {
		assert( that.size_ == size_ ) ;
		//std::copy( that.data_, that.data_ + that.size_, data_ ) ; -> std::copy doesn't seem to work with boost::shared_ptr ??
		assign( *this, that ) ;
		return *this ;
	}

	dense_vector& operator= ( value_type const& value ) {
		assert( size_ == 1 );
		data_[0] = value ;
		return *this;
	}

	dense_vector& operator= ( std::initializer_list< value_type > const& that ) {
		assert( that.size() == size_ );
		//std::copy( that.begin(), that.end(), data_ ) ; -> std::copy doesn't seem to work with boost::shared_ptr ??
		size_type i = 0 ;
		for ( auto x = that.begin(); x != that.end(); ++x, ++i) {
			data_[i] = *x ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, dense_vector& >::type operator= ( E const& that ) {
		assert( that.size() == size_ );
		size_type i = 0 ;
		for ( auto x = that.begin(); x != that.end(); ++x, ++i) {
			data_[i] = *x ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, dense_vector& >::type operator= ( E const& that ) {
		assert( that.size() == size_ ) ;
		assign( *this, that ) ;
		return *this ;
	}

private:
	// Constructor used in shallow_copy
	dense_vector ( boost::shared_ptr<value_type[]> data, boost::shared_ptr<shape_type> shape, size_type size )
	: data_( data )
	, shape_( shape )
	, size_( size )
	{}

private:
	boost::shared_ptr<value_type[]>  data_ ;
	boost::shared_ptr<shape_type>    shape_ ;
	size_type                        size_ ;

public:
	dense_vector shallow_copy () const {
		return dense_vector( data_, shape_, size_ );
	}

public:
	boost::shared_ptr<value_type[]> const data_ptr () const {
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
		return data_[i] ;
	}

	template <typename X>
	auto operator[] ( X const& x ) const -> typename std::enable_if< is<Array, X>::value, decltype( linear_index_selection( *this, x ) ) >::type {
		return linear_index_selection( *this, x ) ;
	}

	template <typename I>
	auto operator[] ( std::initializer_list<I> const& x ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( linear_index_selection( *this, x ) ) >::type {
		return linear_index_selection( *this, x ) ;
	}

	template <typename X>
	typename std::enable_if< is<Vector, X>::value, value_type& >::type operator() ( X const& x ) const {
		assert( x.size() == 1 ) ;
		assert( x[0] >= 0 && x[0] < size_ ) ;
		return data_[x[0]] ;
	}

	value_type&	operator() ( std::initializer_list<size_type> const& x ) const {
		assert( x.size() == 1 ) ;
		assert( *x.begin() >= 0 && *x.begin() < size_ ) ;
		return data_[*x.begin()] ;
	}

	template <typename = void>
	auto operator() ( std::vector<primitive_vector_wrapper<size_type>> const& s ) const -> decltype ( block_selection_of_vector( *this, *s.begin() ) ) {
		assert( s.size() == 1 ) ;
		return block_selection_of_vector( *this, *s.begin() ) ;
	}

} ;

template <typename V>
struct concept< dense_vector<V>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: ContiguousDenseVector
  {};

} ; // namespace glas3

#endif
