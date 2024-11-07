////	  (C) Copyright Sam Corveleyn 2015.
////  Use, modification and distribution are subject to the
////  GLAS Software License, Version 1.0. (See accompanying file
////  LICENSE_1_0.txt)
//
#ifndef glas3_array_dense_array_type_ttt_1D_operation_hpp
#define glas3_array_dense_array_type_ttt_1D_operation_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/type/ttt_1D_operation.hpp>
#include <glas3/array/type/primitive_vector_wrapper.hpp>
#include <glas3/array/dense_array/type/shape_index.hpp>
#include <glas3/array/dense_array/type/range.hpp>

#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection.hpp>

#include <initializer_list>
#include <vector>
#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

namespace glas3 {

template < typename X, typename Y >
class ttt_1D_operation< X, Y, typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value >::type > {

public:
	typedef dense_vector<typename X::size_type>                                                        shape_type ;
	typedef typename shape_type::value_type                                                            size_type ;
	typedef typename shape_type::size_type                                                             ndims_type ;
	typedef decltype( typename X::value_type() * typename Y::value_type() )                            value_type ;

public:
	ttt_1D_operation ( )
	: x_( boost::make_shared<X>( ) )
	, y_( boost::make_shared<Y>( ) )
	, inner_dim_x_( )
	, inner_dim_y_( )
	, outer_dims_x1_( boost::make_shared<shape_type> ( ) )
	, outer_dims_x2_( boost::make_shared<shape_type> ( ) )
	, outer_dims_y_( boost::make_shared<shape_type> ( ) )
	, outer_index_( boost::make_shared<shape_index<typename X::size_type>>( ) )
	, shape_( boost::make_shared<shape_type> ( ) )
	, size_( 0 )
    {}

	ttt_1D_operation ( X const& x, Y const& y, ndims_type const& inner_dim_x, ndims_type const& inner_dim_y )
	: x_( boost::make_shared<X>( std::move( x.shallow_copy() ) ) )
	, y_( boost::make_shared<Y>( std::move( y.shallow_copy() ) ) )
	, inner_dim_x_( inner_dim_x )
	, inner_dim_y_( inner_dim_y )
	, outer_dims_x1_( boost::make_shared<shape_type> ( no_init(), inner_dim_x ) )
	, outer_dims_x2_( boost::make_shared<shape_type> ( no_init(), x.shape().size() - inner_dim_x - 1 ) )
	, outer_dims_y_( boost::make_shared<shape_type> ( no_init(), y.shape().size() - 1 ) )
	, outer_index_( )
	, shape_( boost::make_shared<shape_type> ( no_init(), x.shape().size() + y.shape().size() - 2 ) )
	, size_( 1 )
	{
		assert( x_->shape()[inner_dim_x_] == y_->shape()[inner_dim_y_] ) ;
		size_type i, j = 0, k = 0;

		for ( i = 0; i < inner_dim_x_; ++i, ++j, ++k ) {
				(*shape_)[j] = x_->shape()[i] ;
				size_ *= (*shape_)[j] ;
				(*outer_dims_x1_)[k] = i ;
		}
		k = 0 ;
		for ( i = 0; i < inner_dim_y_; ++i, ++j, ++k ) {
				(*shape_)[j] = y_->shape()[i] ;
				size_ *= (*shape_)[j] ;
				(*outer_dims_y_)[k] = i ;
		}
		for ( i = inner_dim_y_ + 1; i < y_->shape().size(); ++i, ++j, ++k ) {
				(*shape_)[j] = y_->shape()[i] ;
				size_ *= (*shape_)[j] ;
				(*outer_dims_y_)[k] = i ;
		}
		k = 0 ;
		for ( i = inner_dim_x_ + 1; i < x_->shape().size(); ++i, ++j, ++k ) {
				(*shape_)[j] = x_->shape()[i] ;
				size_ *= (*shape_)[j] ;
				(*outer_dims_x2_)[k] = i ;
		}
		outer_index_ = boost::make_shared<shape_index<typename X::size_type>>( *shape_ ) ;

	}

public:
	// Copy constructor -> deep copy
	ttt_1D_operation ( ttt_1D_operation const& that )
	: x_( boost::make_shared<X>( *that.x_ ) )
	, y_( boost::make_shared<Y>( *that.y_ ) )
	, inner_dim_x_( that.inner_dim_x_ )
	, inner_dim_y_( that.inner_dim_y_ )
	, outer_dims_x1_( boost::make_shared<shape_type> ( *that.outer_dims_x1_ ) )
	, outer_dims_x2_( boost::make_shared<shape_type> ( *that.outer_dims_x2_ ) )
	, outer_dims_y_( boost::make_shared<shape_type> ( *that.outer_dims_y_ ) )
	, outer_index_( boost::make_shared<shape_index<typename X::size_type>>( *that.outer_index_ ) )
	, shape_( boost::make_shared<shape_type> ( *that.shape_ ) )
	, size_( that.size_ )
	{}

	// Move constructor
	ttt_1D_operation ( ttt_1D_operation && that ) = default ;

public:
	// Copy assignment -> do not allow
	ttt_1D_operation& operator= ( ttt_1D_operation const& that ) = delete ; // not assignable

	// Move assignment
	ttt_1D_operation& operator= ( ttt_1D_operation&& that ) = delete ;

private:
	// Constructor used in shallow_copy
	ttt_1D_operation ( boost::shared_ptr< X > x, boost::shared_ptr< Y > y,	ndims_type inner_dim_x,	ndims_type inner_dim_y,
			boost::shared_ptr< shape_type > outer_dims_x1, boost::shared_ptr< shape_type > outer_dims_x2, boost::shared_ptr< shape_type > outer_dims_y,
			boost::shared_ptr< shape_index<typename X::size_type> > outer_index,
			boost::shared_ptr< shape_type > shape, size_type size )
	: x_( x )
	, y_( y )
	, inner_dim_x_( inner_dim_x )
	, inner_dim_y_( inner_dim_y )
	, outer_dims_x1_( outer_dims_x1 )
	, outer_dims_x2_( outer_dims_x2 )
	, outer_dims_y_( outer_dims_y )
	, outer_index_( outer_index )
	, shape_( shape )
	, size_( size )
	{}

private:
	boost::shared_ptr< X >                                                       x_ ;
	boost::shared_ptr< Y >                                                       y_ ;
	ndims_type                                                                   inner_dim_x_ ;
	ndims_type                                                                   inner_dim_y_ ;
	boost::shared_ptr< shape_type >                                              outer_dims_x1_ ;
	boost::shared_ptr< shape_type >                                              outer_dims_x2_ ;
	boost::shared_ptr< shape_type >                                              outer_dims_y_ ;
	boost::shared_ptr< shape_index<typename X::size_type> >                      outer_index_ ;
	boost::shared_ptr< shape_type >                                              shape_ ;
	size_type                                                                    size_ ;

public:
	ttt_1D_operation shallow_copy () const {
		return ttt_1D_operation( x_, y_, inner_dim_x_, inner_dim_y_, outer_dims_x1_, outer_dims_x2_, outer_dims_y_, outer_index_, shape_, size_ );
	}

public:
	shape_type const& shape () const {
		return *shape_ ;
	}

	size_type size () const {
		return size_ ;
	}

	size_type ndof () const {
		return x_->ndof() + y_->ndof() ;
	}

public:
    template <typename I>
	typename std::enable_if< std::is_integral<I>::value, value_type >::type operator[] ( I const& i ) const {
		outer_index_->set_lin_in( i ) ;
		dense_vector<typename X::size_type> j_x( no_init(), x_->shape().size() ) ;
		dense_vector<typename Y::size_type> j_y( no_init(), y_->shape().size() ) ;
		j_x[*outer_dims_x1_] = (*outer_index_)[range<>( 0, inner_dim_x_ )] ;
		j_y[*outer_dims_y_] = (*outer_index_)[range<>( inner_dim_x_, inner_dim_x_ + y_->shape().size() - 1 )] ;
		j_x[*outer_dims_x2_] = (*outer_index_)[range<>( inner_dim_x_ + y_->shape().size() - 1, shape_->size() )] ;
		value_type sum = 0 ;
		for ( size_type j = 0; j < x_->shape()[inner_dim_x_]; ++j ) {
			j_x[inner_dim_x_] = j ;
			j_y[inner_dim_y_] = j ;
			sum += (*x_)(j_x) * (*y_)(j_y) ;
		}
		return sum ;
	}

    template <typename S>
    auto operator[] ( S const& x ) const -> typename std::enable_if< is<Array, S>::value, decltype ( linear_index_selection( *this, x ) ) >::type {
    	return linear_index_selection( *this, x ) ;
    }

    template <typename I>
    auto operator[] ( std::initializer_list<I> const& x ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( linear_index_selection( *this, x ) ) >::type {
      	return linear_index_selection( *this, x ) ;
    }

	template <typename S>
	typename std::enable_if< is<Vector, S>::value, value_type >::type operator() ( S const& j ) const {
		assert( j.size() == shape_->size() ) ;
		dense_vector<typename X::size_type> j_x( no_init(), x_->shape().size() ) ;
		dense_vector<typename Y::size_type> j_y( no_init(), y_->shape().size() ) ;
		j_x[*outer_dims_x1_] = j[range<>( 0, inner_dim_x_ )] ;
		j_y[*outer_dims_y_] = j[range<>( inner_dim_x_, inner_dim_x_ + y_->shape().size() - 1 )] ;
		j_x[*outer_dims_x2_] = j[range<>( inner_dim_x_ + y_->shape().size() - 1, shape_->size() )] ;
		value_type sum = 0 ;
		for ( size_type j = 0; j < x_->shape()[inner_dim_x_]; ++j ) {
			j_x[inner_dim_x_] = j ;
			j_y[inner_dim_y_] = j ;
			sum += (*x_)(j_x) * (*y_)(j_y) ;
		}
		return sum ;
	}

	value_type operator() ( std::initializer_list<size_type> const& x ) const {
		return (*this)( dense_vector<size_type>( x ) ) ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

} ;

template < typename X, typename Y >
struct concept< ttt_1D_operation< X, Y >, typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value
                                                            && ! ( is< DenseVector, X >::value && is< DenseVector, Y >::value )
                                                            && ! ( is< DenseMatrix, X >::value && is< DenseVector, Y >::value )
                                                            && ! ( is< DenseVector, X >::value && is< DenseMatrix, Y >::value )
                                                            && ! ( is< DenseMatrix, X >::value && is< DenseMatrix, Y >::value ) >::type >
: DenseArray
  {} ;

template < typename X, typename Y >
struct concept< ttt_1D_operation< X, Y >, typename std::enable_if< is< DenseVector, X >::value && is< DenseVector, Y >::value >::type >
: DenseScalar
  {} ;

template < typename X, typename Y >
struct concept< ttt_1D_operation< X, Y >, typename std::enable_if< is< DenseMatrix, X >::value && is< DenseVector, Y >::value >::type >
: DenseVector
  {} ;

template < typename X, typename Y >
struct concept< ttt_1D_operation< X, Y >, typename std::enable_if< is< DenseVector, X >::value && is< DenseMatrix, Y >::value >::type >
: DenseVector
  {} ;

template < typename X, typename Y >
struct concept< ttt_1D_operation< X, Y >, typename std::enable_if< is< DenseMatrix, X >::value && is< DenseMatrix, Y >::value >::type >
: DenseMatrix
  {} ;

} ; // namespace glas3

#endif
