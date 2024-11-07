//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_eye_hpp
#define glas3_array_dense_array_type_eye_hpp

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
class eye {
	static_assert( std::is_arithmetic<V>::value || boost::is_complex<V>::value, "V should be arithmetic" ) ;
public:
    typedef V                                 value_type ;
    typedef dense_vector<std::ptrdiff_t>      shape_type ;
    typedef typename shape_type::value_type   size_type ;
    typedef typename shape_type::size_type    ndims_type ;

public:
    eye( )
    : shape_( boost::make_shared< shape_type >( ) )
	, dim_factor_( boost::make_shared< shape_type > ( ) )
	, j_mem_ ( boost::make_shared< shape_type >( ) )
	, i_in_ ( boost::make_shared< size_type >( ) )
	, i_out_ ( boost::make_shared< size_type >( ) )
    , size_( 1 )
    {
    	if ( ndims >= 0 ) {
    		assert ( ndims == 0 ) ;
    	}
    }

    eye( size_type dim_size, ndims_type d )
    : shape_( boost::make_shared< shape_type >( dim_size, d ) )
	, dim_factor_( boost::make_shared< shape_type > ( 1, d ) )
	, j_mem_ ( boost::make_shared< shape_type >( 0, d ) )
	, i_in_ ( boost::make_shared< size_type >( 0 ) )
	, i_out_ ( boost::make_shared< size_type >( 0 ) )
    , size_( std::pow( dim_size, d ) )
    {
    	if ( ndims >= 0 ) {
    		assert ( ndims == d ) ;
    	}
		for ( ndims_type k = 1; k < shape_->size(); ++k  ) {
			(*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ;
		}
    }

    eye( std::initializer_list< size_type > const& shape )
    : shape_( boost::make_shared< shape_type >( std::move ( shape_type ( shape ) ) ) )
	, dim_factor_( boost::make_shared< shape_type > ( 1, shape.size() ) )
	, j_mem_ ( boost::make_shared< shape_type >( 0, shape.size() ) )
	, i_in_ ( boost::make_shared< size_type >( 0 ) )
	, i_out_ ( boost::make_shared< size_type >( 0 ) )
    , size_( 1 )
    {
    	if ( ndims >= 0 ) {
			assert ( ndims == shape_->size() ) ;
		}
		for ( ndims_type k = 0; k < shape_->size(); ++k ){
			size_ *= (*shape_)[k] ;
			if ( k > 0 ) { (*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ; }
		}
    }

    eye( shape_type shape )
    : shape_( boost::make_shared< shape_type >( std::move ( shape ) ) )
	, dim_factor_( boost::make_shared< shape_type > ( 1, shape_->size() ) )
	, j_mem_ ( boost::make_shared< shape_type >( 0, shape_->size() ) )
	, i_in_ ( boost::make_shared< size_type >( 0 ) )
	, i_out_ ( boost::make_shared< size_type >( 0 ) )
    , size_( 1 )
    {
    	if ( ndims >= 0 ) {
    		assert ( ndims == shape_->size() ) ;
    	}
    	for ( ndims_type k = 0; k < shape_->size(); ++k ){
			size_ *= (*shape_)[k] ;
			if ( k > 0 ) { (*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ; }
		}
    }

public:
    // Copy constructor -> deep copy
    eye ( eye const& that )
	: shape_( boost::make_shared< shape_type >( *that.shape_) )
    , dim_factor_( boost::make_shared< shape_type >( *that.dim_factor_ ) )
	, j_mem_ ( boost::make_shared< shape_type >( *that.j_mem_ ) )
	, i_in_ ( boost::make_shared< size_type >( *that.i_in_ ) )
	, i_out_ ( boost::make_shared< size_type >( *that.i_out_ ) )
	, size_( that.size_ )
    {}

    // Move constructor
    eye ( eye && that ) = default ;

public:
    // Copy assignment -> do not allow
    eye& operator= ( eye const& that ) = delete ;

    // Move assignment
    eye& operator= ( eye&& that ) = delete ;

private:
    // Constructor used in shallow_copy
	eye ( boost::shared_ptr<shape_type> shape, boost::shared_ptr< shape_type >  dim_factor, boost::shared_ptr< shape_type > j_mem,
			boost::shared_ptr< size_type > i_in, boost::shared_ptr< size_type > i_out, size_type size )
    : shape_( shape )
	, dim_factor_( dim_factor )
    , j_mem_ ( j_mem )
    , i_in_ ( i_in )
    , i_out_ ( i_out )
    , size_( size )
    {}

private:
    boost::shared_ptr<shape_type>    shape_ ;
	boost::shared_ptr< shape_type >  dim_factor_ ;
	boost::shared_ptr< shape_type >  j_mem_ ;
	boost::shared_ptr< size_type >   i_in_ ;
	boost::shared_ptr< size_type >   i_out_ ;
    size_type                        size_ ;

public:
    eye shallow_copy () const {
       	return eye( shape_, dim_factor_, j_mem_, i_in_, i_out_, size_ );
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
		if ( shape_->size() == 0 ) {
			return 1 ;
		}
		else if ( shape_->size() == 1 ) {
			if ( i > 0 ) { return 0 ; }
			else { return 1 ; }
		}
		else {
			index2multi_index_memory( *dim_factor_, *shape_, *j_mem_, i, *i_in_) ;

			auto j_val0 = (*j_mem_)[0] ;
			for ( ndims_type k = 1; k < shape_->size(); ++k ) {
				if ( (*j_mem_)[k] != j_val0 ) return 0 ;
			}
			return 1 ;
		}
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
    	if ( shape_->size() == 0 ) {
			return 1 ;
		}
		else if ( shape_->size() == 1 ) {
			auto x_val0 = x[0] ;
			assert( x_val0 >= 0 && x_val0 < (*shape_)[0] ) ;
			if ( x_val0 > 0 ) { return 0 ; }
			else { return 1 ; }
		}
		else {
			auto x_val0 = x[0] ;
			assert( x_val0 >= 0 && x_val0 < (*shape_)[0] ) ;
			for ( ndims_type k = 1; k < shape_->size(); ++k ) {
				assert( x[k] >= 0 && x[k] < (*shape_)[k] ) ;
				if ( x[k] != x_val0 ) return 0 ;
			}
			return 1 ;
		}
    }

	value_type operator() ( std::initializer_list<size_type> const& x ) const {
    	assert( x.size() == shape_->size() ) ;
    	if ( shape_->size() == 0 ) {
    		return 1 ;
    	}
    	else if ( shape_->size() == 1 ) {
    		auto x_val0 = *x.begin() ;
    		assert( x_val0 >= 0 && x_val0 < (*shape_)[0] ) ;
    		if ( x_val0 > 0 ) { return 0 ; }
    		else { return 1 ; }
    	}
    	else {
			auto x_it = x.begin() ;
			auto x_val0 = *x_it ; ++x_it ;
			assert( x_val0 >= 0 && x_val0 < (*shape_)[0] ) ;
			ndims_type k = 1 ;
			for ( ; x_it != x.end(); ++x_it, ++k ) {
				assert( *x_it >= 0 && *x_it < (*shape_)[k] ) ;
				if ( *x_it != x_val0 ) return 0 ;
			}
			return 1 ;
    	}
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

  } ;

template <typename V, std::ptrdiff_t ndims>
struct concept< eye<V, ndims>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseArray
{} ;

template <typename V>
struct concept< eye<V, 0>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseScalar
{} ;

template <typename V>
struct concept< eye<V, 1>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseVector
{} ;

template <typename V>
struct concept< eye<V, 2>, typename std::enable_if< std::is_arithmetic<V>::value || boost::is_complex<V>::value >::type >
: DenseMatrix
{} ;

} // namespace glas3

#endif
