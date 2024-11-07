////	  (C) Copyright Sam Corveleyn 2015.
////  Use, modification and distribution are subject to the
////  GLAS Software License, Version 1.0. (See accompanying file
////  LICENSE_1_0.txt)
//
#ifndef glas3_array_dense_array_type_ttt_operation_hpp
#define glas3_array_dense_array_type_ttt_operation_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/type/ttt_operation.hpp>
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

template < typename X, typename Y, typename Dx, typename Dy >
class ttt_operation< X, Y, Dx, Dy, typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value && is< Array, Dx >::value && is< Array, Dy >::value >::type > {

public:
	typedef dense_vector<typename X::size_type>                                                                                 shape_type ;
	typedef typename shape_type::value_type                                                                                     size_type ;
	typedef typename shape_type::size_type                                                                                      ndims_type ;
	typedef decltype( typename X::value_type() * typename Y::value_type() )                                                     value_type ;

public:
	ttt_operation ( )
	: x_( boost::make_shared<X>( ) )
	, y_( boost::make_shared<Y>( ) )
	, inner_dims_x_( boost::make_shared<Dx> ( ) )
	, inner_dims_y_( boost::make_shared<Dy> ( ) )
	, outer_dims_x_( boost::make_shared<shape_type> ( ) )
	, outer_dims_y_( boost::make_shared<shape_type> ( ) )
    , inner_index_( boost::make_shared<shape_index<typename X::size_type>>( ) )
	, outer_index_( boost::make_shared<shape_index<typename X::size_type>>( ) )
	, shape_( boost::make_shared<shape_type> ( ) )
	, size_( 0 )
    {}

	ttt_operation ( X const& x, Y const& y, Dx inner_dims_x, Dy inner_dims_y )
	: x_( boost::make_shared<X>( std::move( x.shallow_copy() ) ) )
	, y_( boost::make_shared<Y>( std::move( y.shallow_copy() ) ) )
	, inner_dims_x_( boost::make_shared<Dx> ( std::move( inner_dims_x ) ) )
	, inner_dims_y_( boost::make_shared<Dy> ( std::move( inner_dims_y ) ) )
	, outer_dims_x_( boost::make_shared<shape_type> ( no_init(), x.shape().size() - inner_dims_x.size() ) )
	, outer_dims_y_( boost::make_shared<shape_type> ( no_init(), y.shape().size() - inner_dims_y.size() ) )
    , inner_index_( boost::make_shared<shape_index<typename X::size_type>>( x_->shape()[*inner_dims_x_] ) )
	, outer_index_( )
	, shape_( boost::make_shared<shape_type> ( no_init(), x.shape().size() + y.shape().size() - inner_dims_x.size() - inner_dims_y.size() ) )
	, size_( 1 )
	{
		assert( inner_dims_x_->size() == inner_dims_y_->size() );
		size_type i, k = 0, l = 0 ;
		shape_type mask_inner_dims_x( 0, x_->shape().size() ), mask_inner_dims_y( 0, y_->shape().size() ) ;

		for ( i = 0; i < inner_dims_x_->size(); ++i ) {
			assert( x_->shape()[(*inner_dims_x_)[i]] == y_->shape()[(*inner_dims_y_)[i]] ) ;
			mask_inner_dims_x[(*inner_dims_x_)[i]] = 1 ;
			mask_inner_dims_y[(*inner_dims_y_)[i]] = 1 ;
		}
		for ( i = 0; i < x_->shape().size(); ++i  ) {
			if ( ! mask_inner_dims_x[i] ) {
				(*shape_)[k] = x_->shape()[i] ;
				size_ *= (*shape_)[k] ;
				(*outer_dims_x_)[k] = i ;
				++k ;
			}
		}
		for ( i = 0; i < y_->shape().size(); ++i ) {
			if ( ! mask_inner_dims_y[i] ) {
				(*shape_)[k] = y_->shape()[i] ;
				size_ *= (*shape_)[k] ;
				(*outer_dims_y_)[l] = i ;
				++k ;
				++l ;
			}
		}
		outer_index_ = boost::make_shared<shape_index<typename X::size_type>>( *shape_ ) ;
	}

public:
	// Copy constructor -> deep copy
	ttt_operation ( ttt_operation const& that )
	: x_( boost::make_shared<X>( *that.x_ ) )
	, y_( boost::make_shared<Y>( *that.y_ ) )
	, inner_dims_x_( boost::make_shared<Dx> ( *that.inner_dims_x_ ) )
	, inner_dims_y_( boost::make_shared<Dy> ( *that.inner_dims_y_ ) )
	, outer_dims_x_( boost::make_shared<shape_type> ( *that.outer_dims_x_ ) )
	, outer_dims_y_( boost::make_shared<shape_type> ( *that.outer_dims_y_ ) )
	, inner_index_( boost::make_shared<shape_index<typename X::size_type>>( *that.inner_index_ ) )
	, outer_index_( boost::make_shared<shape_index<typename X::size_type>>( *that.outer_index_ ) )
	, shape_( boost::make_shared<shape_type> ( *that.shape_ ) )
	, size_( that.size_ )
	{}

	// Move constructor
	ttt_operation ( ttt_operation && that ) = default;

public:
	// Copy assignment -> do not allow
	ttt_operation& operator= ( ttt_operation const& that ) = delete ; // not assignable

	// Move assignment
	ttt_operation& operator= ( ttt_operation&& that ) = delete ;

private:
	// Constructor used in shallow_copy
	ttt_operation ( boost::shared_ptr< X > x, boost::shared_ptr< Y > y,
			boost::shared_ptr< Dx > inner_dims_x, boost::shared_ptr< Dy > inner_dims_y,
			boost::shared_ptr< shape_type > outer_dims_x, boost::shared_ptr< shape_type > outer_dims_y,
			boost::shared_ptr< shape_index<typename X::size_type> > inner_index, boost::shared_ptr< shape_index<typename X::size_type> > outer_index,
			boost::shared_ptr< shape_type > shape, size_type size )
	: x_( x )
	, y_( y )
	, inner_dims_x_( inner_dims_x )
	, inner_dims_y_( inner_dims_y )
	, outer_dims_x_( outer_dims_x )
	, outer_dims_y_( outer_dims_y )
	, inner_index_( inner_index )
	, outer_index_( outer_index )
	, shape_( shape )
	, size_( size )
	{}

private:
	boost::shared_ptr< X >                                                       x_ ;
	boost::shared_ptr< Y >                                                       y_ ;
	boost::shared_ptr< Dx >                                                      inner_dims_x_ ;
	boost::shared_ptr< Dy >                                                      inner_dims_y_ ;
	boost::shared_ptr< shape_type >                                              outer_dims_x_ ;
	boost::shared_ptr< shape_type >                                              outer_dims_y_ ;
	boost::shared_ptr< shape_index<typename X::size_type> >                      inner_index_ ;
	boost::shared_ptr< shape_index<typename X::size_type> >                      outer_index_ ;
	boost::shared_ptr< shape_type >                                              shape_ ;
	size_type                                                                    size_ ;

public:
	ttt_operation shallow_copy () const {
		return ttt_operation( x_, y_, inner_dims_x_, inner_dims_y_, outer_dims_x_, outer_dims_y_, inner_index_, outer_index_, shape_, size_ );
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
		j_x[*outer_dims_x_] = (*outer_index_)[range<>( 0, outer_dims_x_->size() )] ;
		j_y[*outer_dims_y_] = (*outer_index_)[range<>( outer_dims_x_->size(), shape_->size() )] ;
		value_type sum = 0 ;
		for ( inner_index_->reset(); inner_index_->overflow_count() < 1; inner_index_->inc_lin_in() ) {
			j_x[*inner_dims_x_] = *inner_index_ ;
			j_y[*inner_dims_y_] = *inner_index_ ;
			sum += (*x_)(j_x) * (*y_)(j_y) ;
		}
		return sum ;
	}

	template <typename S>
	auto operator[] ( S const& s ) const -> typename std::enable_if< is<Array, S>::value, decltype( linear_index_selection( *this, s ) ) >::type {
		return linear_index_selection( *this, s ) ;
	}

    template <typename I>
    auto operator[] ( std::initializer_list<I> const& x ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( linear_index_selection( *this, x ) ) >::type {
      	return linear_index_selection( *this, x ) ;
    }

	template <typename S>
	typename std::enable_if< is<Vector, S>::value, value_type >::type
	operator() ( S const& j ) const {
		assert( j.size() == shape_->size() ) ;
		dense_vector<typename X::size_type> j_x( no_init(), x_->shape().size() ) ;
		dense_vector<typename Y::size_type> j_y( no_init(), y_->shape().size() ) ;
		j_x[*outer_dims_x_] = j[range<>( 0, outer_dims_x_->size() )] ;
		j_y[*outer_dims_y_] = j[range<>( outer_dims_x_->size(), shape_->size() )] ;
		value_type sum = 0 ;
		for ( inner_index_->reset(); inner_index_->overflow_count() < 1; inner_index_->inc_lin_in() ) {
			j_x[*inner_dims_x_] = *inner_index_ ;
			j_y[*inner_dims_y_] = *inner_index_ ;
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

template < typename X, typename Y, typename Dx, typename Dy >
struct concept< ttt_operation< X, Y, Dx, Dy >, typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value && is< Array, Dx >::value && is< Array, Dy >::value >::type >
: DenseArray
  {} ;

} ; // namespace glas3

#endif
