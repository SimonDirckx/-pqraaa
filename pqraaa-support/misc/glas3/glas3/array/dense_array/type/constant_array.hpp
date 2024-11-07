//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_constant_array_hpp
#define glas3_array_dense_array_type_constant_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/type/primitive_vector_wrapper.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection.hpp>
#include <glas3/array/algorithm/indexing.hpp>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <vector>
#include <initializer_list>
#include <utility>

namespace glas3 {

template < typename V, std::ptrdiff_t ndims = -1 >
class constant_array {
	static_assert( std::is_arithmetic<V>::value || boost::is_complex<V>::value, "V should be arithmetic" ) ;
public:
    typedef V                                 value_type ;
    typedef dense_vector<std::ptrdiff_t>      shape_type ;
    typedef typename shape_type::value_type   size_type ;
    typedef typename shape_type::size_type    ndims_type ;

public:
    constant_array( )
    : shape_( boost::make_shared< shape_type >( ) )
    , value_( )
    , size_( 0 )
    {
    	if ( ndims >= 0 ) {
    		assert ( ndims == 0 ) ;
    	}
    }

    constant_array ( value_type const& value )
    : shape_( boost::make_shared< shape_type >( ) )
    , value_( value )
    , size_( 1 )
    {
    	if ( ndims >= 0 ) {
    		assert ( ndims == 0 ) ;
    	}
    }

    constant_array( value_type const& value, shape_type shape )
    : shape_( boost::make_shared< shape_type >( std::move ( shape ) ) )
    , value_( value )
    , size_( 1 )
    {
    	if ( ndims >= 0 ) {
    		assert ( ndims == shape_->size() ) ;
    	}
    	for ( ndims_type k = 0; k < shape_->size(); ++k ){
    		size_ *= (*shape_)[k] ;
    	}
    }

public:
    // Copy constructor -> deep copy
    constant_array ( constant_array const& that )
	: shape_( boost::make_shared< shape_type >( *that.shape_ ) )
    , value_( that.value_ )
	, size_( that.size_ )
    {}

    // Move constructor
    constant_array ( constant_array && that ) = default ;

public:
    // Copy assignment -> do not allow
    constant_array& operator= ( constant_array const& that ) = delete ;

    // Move assignment
    constant_array& operator= ( constant_array&& that ) = delete ;

private:
    // Constructor used in shallow_copy
	constant_array ( boost::shared_ptr<shape_type> shape, value_type value, size_type size )
    : shape_( shape )
	, value_( value )
    , size_( size )
    {}

private:
    boost::shared_ptr<shape_type>    shape_ ;
    value_type                       value_ ;
    size_type                        size_ ;

public:
    constant_array shallow_copy () const {
       	return constant_array( shape_, value_, size_ );
    }

public:
    shape_type const& shape () const {
        return *shape_ ;
    }

    size_type size () const {
        return size_ ;
    }

    size_type ndof () const {
        return 0 ;
    }

public:
        template <typename I>
	typename std::enable_if< std::is_integral<I>::value, value_type >::type operator[] ( I const& i ) const {
		assert( i >= 0 && i < size_ ) ;
		return value_ ;
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
    typename std::enable_if< is<Vector, X>::value, value_type >::type operator() ( X const& x ) const {
    	assert( x.size() == shape_->size() ) ;
    	for ( ndims_type k = 0; k < shape_->size(); ++k ) {
    		assert( x[k] >= 0 && x[k] < (*shape_)[k] ) ;
    	}
    	return value_ ;
    }

	value_type operator() ( std::initializer_list<size_type> const& x ) const {
    	assert( x.size() == shape_->size() ) ;
    	ndims_type k = 0 ;
    	for ( auto x_it = x.begin(); x_it != x.end(); ++x_it, ++k ) {
    		assert( *x_it >= 0 && *x_it < (*shape_)[k] ) ;
    	}
    	return value_ ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

  } ;

template <typename V, std::ptrdiff_t ndims>
struct concept< constant_array<V, ndims>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseArray
{} ;

template <typename V>
struct concept< constant_array<V, 0>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseScalar
{} ;

template <typename V>
struct concept< constant_array<V, 1>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseVector
{} ;

template <typename V>
struct concept< constant_array<V, 2>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseMatrix
{} ;

} // namespace glas3

#endif
