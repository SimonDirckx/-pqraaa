//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_range_hpp
#define glas3_array_dense_array_type_range_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/type/primitive_vector_wrapper.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection_of_vector.hpp>

#include <initializer_list>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <vector>
#include <algorithm>

namespace glas3 {

template < typename V = std::ptrdiff_t >
class range {
public:
    typedef V                                 value_type ;
    typedef dense_scalar<std::ptrdiff_t>      shape_type ;
    typedef typename shape_type::value_type   size_type ;
    typedef typename shape_type::size_type    ndims_type ;

public:
    range( )
    : begin_( 0 )
    , end_( 0 )
    , shape_( boost::make_shared< shape_type >( 0 ) )
    , size_( 0 )
    , step_( 1 )
    {}

    range( value_type begin, value_type end )
    : begin_( begin )
    , end_( end )
    , shape_( boost::make_shared< shape_type >( std::max( V(0), end - begin ) ) )
    , size_( std::max( V(0), end - begin ) )
    , step_( 1 )
    {}

    range( value_type begin, value_type end, value_type step )
    : begin_( begin )
    , end_( end )
    , shape_( boost::make_shared< shape_type >( std::max( 0.0, ceil( double( end - begin ) / step ) ) ) )
    , size_( std::max( 0.0, ceil( double( end - begin ) / step ) ) )
    , step_( step )
    {}

public:
    // Copy constructor -> deep copy
    range ( range const& that )
    : begin_( that.begin_ )
    , end_( that.end_ )
    , shape_( boost::make_shared< shape_type >( *that.shape_ ) )
    , size_( that.size_ )
    , step_( that.step_ )
    {}

    // Move constructor
    range ( range && that ) = default ;

public:
    // Copy assignment -> deep copy
    range& operator= ( range const& that ) = delete ;

    // Move assignment
    range& operator= ( range&& that ) = delete ;

private:
    // Constructor used in shallow_copy
	range ( value_type begin, value_type end, boost::shared_ptr<shape_type> shape, size_type size, value_type step )
    : begin_( begin )
    , end_( end )
	, shape_( shape )
    , size_( size )
	, step_( step )
    {}

private:
    value_type                       begin_ ;
    value_type                       end_ ;
    boost::shared_ptr<shape_type>    shape_ ;
    size_type                        size_ ;
    value_type                       step_ ;

public:
    range shallow_copy () const {
       	return range( begin_, end_, shape_, size_, step_ );
    }

public:
    shape_type const& shape () const {
        return *shape_ ;
    }

    size_type size () const {
        return size_ ;
    }

    size_type ndof () const {
        return 3 ;
    }

public:
    value_type begin() const {
    	return begin_ ;
    }

    value_type end() const {
    	return end_ ;
    }

    value_type step() const {
    	return step_ ;
    }

public:
    template <typename I>
    typename std::enable_if< std::is_integral<I>::value, value_type >::type operator[] ( I const& i ) const {
		assert( i >= 0 && i < size_ ) ;
		return begin_ + i * step_ ;
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
    	assert( x.size() == 1 ) ;
    	assert( x[0] >= 0 && x[0] < size_ ) ;
    	return begin_ + x[0] * step_ ;
    }

	value_type operator() ( std::initializer_list<size_type> const& x ) const {
    	assert( x.size() == 1 ) ;
    	assert( *x.begin() >= 0 && *x.begin() < size_ ) ;
    	return begin_ + (*x.begin()) * step_ ;
	}

	template <typename = void>
	auto operator() ( std::vector<primitive_vector_wrapper<size_type>> const& s ) const -> decltype ( block_selection_of_vector( *this, *s.begin() ) ) {
		assert( s.size() == 1 ) ;
		return block_selection_of_vector( *this, *s.begin() ) ;
	}

  } ;

template <typename V>
struct concept< range<V>, typename std::enable_if< std::is_integral<V>::value >::type >
: DenseVector
{} ;

} // namespace glas3

#endif
