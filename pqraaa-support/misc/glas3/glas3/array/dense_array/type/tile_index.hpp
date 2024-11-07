//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_tile_index_hpp
#define glas3_array_dense_array_type_tile_index_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/dense_array/concept/index.hpp>

#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/type/primitive_vector_wrapper.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection_of_vector.hpp>
#include <glas3/array/algorithm/indexing.hpp>

#include <initializer_list>
#include <type_traits>
#include <cassert>
#include <vector>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

namespace glas3 {

template <typename V = std::ptrdiff_t>
class tile_index {

public:
	typedef V                                  value_type;
	typedef dense_vector<value_type>           index_type ;
	typedef dense_scalar<std::ptrdiff_t>       shape_type ;
	typedef typename shape_type::value_type    size_type ;
	typedef typename shape_type::size_type     ndims_type ;

public:
	tile_index ( )
	: shape_in_( boost::make_shared< index_type >( ) )
    , shape_out_( boost::make_shared< index_type >( ) )
	, tile_shape_( boost::make_shared< index_type >( ) )
	, dim_factor_in_( boost::make_shared< index_type >( ) )
	, dim_factor_out_( boost::make_shared< index_type >( ) )
	, j_in_( boost::make_shared< index_type >( ) )
	, j_out_( boost::make_shared< index_type >( ) )
	, lin_size_( )
    , overflow_count_( )
	, i_in_( boost::make_shared< size_type >( ) )
	, i_out_( boost::make_shared< size_type >( ) )
	, i_out_of_date_( false )
	, shape_( boost::make_shared< shape_type >( 0 ) )
	, size_( 0 )
	{}

	tile_index ( index_type shape_out, index_type tile_shape )
	: shape_in_( boost::make_shared< index_type >( no_init(), shape_out.size() ) )
    , shape_out_( boost::make_shared< index_type >( std::move( shape_out ) ) )
	, tile_shape_( boost::make_shared< index_type >( std::move( tile_shape ) ) )
	, dim_factor_in_( boost::make_shared< index_type >( 1, shape_in_->size() ) )
	, dim_factor_out_( boost::make_shared< index_type >( 1, shape_out_->size() ) )
	, j_in_( boost::make_shared< index_type >( 0, shape_in_->size() ) )
	, j_out_( boost::make_shared< index_type >( 0, shape_out_->size() ) )
	, lin_size_( 1 )
    , overflow_count_( 0 )
	, i_in_( boost::make_shared< size_type >( 0 ) )
	, i_out_( boost::make_shared< size_type >( 0 ) )
	, i_out_of_date_( false )
	, shape_( boost::make_shared< shape_type >( shape_out_->size() ) )
	, size_( shape_out_->size() )
	{
	    for ( size_type k = 0; k < shape_out_->size(); ++k ) {
	    	(*shape_in_)[k] = (*shape_out_)[k] * (*tile_shape_)[k] ;
	    	lin_size_ *= (*shape_in_)[k] ;
		    if ( k > 0 ) {
			    (*dim_factor_in_)[k] = (*dim_factor_in_)[k - 1] * (*shape_in_)[k - 1];
			    (*dim_factor_out_)[k] = (*dim_factor_out_)[k - 1] * (*shape_out_)[k - 1];
		    }
	    }
	}

public:
	// Copy constructor -> deep copy
	tile_index ( tile_index const& that )
	: shape_in_( boost::make_shared< index_type >( *that.shape_in_ ) )
	, shape_out_( boost::make_shared< index_type >( *that.shape_out_ ) )
	, tile_shape_( boost::make_shared< index_type >( *that.tile_shape_ ) )
	, dim_factor_in_( boost::make_shared< index_type >( *that.dim_factor_in_ ) )
	, dim_factor_out_( boost::make_shared< index_type >( *that.dim_factor_out_ ) )
	, j_in_( boost::make_shared< index_type >( *that.j_in_ ) )
	, j_out_( boost::make_shared< index_type >( *that.j_out_ ) )
	, lin_size_( that.lin_size_ )
	, overflow_count_( that.overflow_count_ )
	, i_in_( boost::make_shared< size_type >( *that.i_in_ ) )
	, i_out_( boost::make_shared< size_type >( *that.i_out_ ) )
	, i_out_of_date_( that.i_out_of_date_ )
	, shape_( boost::make_shared< shape_type >( *that.shape_ ) )
	, size_( that.size_ )
	{}

	// Move constructor
	tile_index ( tile_index && that ) = default ;

public:
	// Copy assignment -> do not allow
	tile_index& operator= ( tile_index const& that ) = delete ; // not assignable

	// Move assignment
	tile_index& operator= ( tile_index&& that) = delete ;

private:
	// Constructor used in shallow_copy
	tile_index ( boost::shared_ptr<index_type> shape_in, boost::shared_ptr<index_type> shape_out, boost::shared_ptr<index_type> tile_shape,
			boost::shared_ptr<index_type> dim_factor_in, boost::shared_ptr<index_type> dim_factor_out,
			boost::shared_ptr<index_type> j_in, boost::shared_ptr<index_type> j_out,
			value_type lin_size, value_type overflow_count,
			boost::shared_ptr<value_type> i_in, boost::shared_ptr<value_type> i_out, bool i_out_of_date,
			boost::shared_ptr<shape_type> shape, size_type size)
	: shape_in_( shape_in )
	, shape_out_( shape_out )
	, tile_shape_( tile_shape )
	, dim_factor_in_( dim_factor_in )
	, dim_factor_out_( dim_factor_out )
	, j_in_( j_in )
	, j_out_( j_out )
	, lin_size_( lin_size )
	, overflow_count_( overflow_count )
	, i_in_( i_in )
	, i_out_( i_out )
	, i_out_of_date_( i_out_of_date )
	, shape_( shape )
	, size_( size )
	{}

private:
	boost::shared_ptr<index_type>      shape_in_ ;
	boost::shared_ptr<index_type>      shape_out_ ;
	boost::shared_ptr<index_type>      tile_shape_ ;
	boost::shared_ptr<index_type>      dim_factor_in_ ;
	boost::shared_ptr<index_type>      dim_factor_out_ ;
	boost::shared_ptr<index_type>      j_in_ ;
	boost::shared_ptr<index_type>      j_out_ ;
	value_type                         lin_size_ ;
	mutable value_type                 overflow_count_ ;
	boost::shared_ptr<value_type>      i_in_ ;
	boost::shared_ptr<value_type>      i_out_ ;
	mutable bool                       i_out_of_date_ ;
	boost::shared_ptr<shape_type>      shape_ ;
	size_type                          size_ ;

public:
	tile_index shallow_copy () const {
		return tile_index( shape_in_, shape_out_, tile_shape_, dim_factor_in_, dim_factor_out_,
				j_in_, j_out_, lin_size_, overflow_count_, i_in_, i_out_, i_out_of_date_, shape_, size_ );
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

private:
	void update_i () const {
		*i_in_ = multi_index2index_initializer_list( *shape_in_, *j_in_ ) ;
		*i_out_ = multi_index2index_initializer_list( *shape_out_, *j_out_ ) ;
		i_out_of_date_ = false ;
	}

public:
	value_type lin_size () const {
		return lin_size_ ;
	}

	value_type overflow_count () const {
		return overflow_count_ ;
	}

	void set_lin_in ( value_type const& i ) const {
		assert( i >= 0 && i < lin_size_ ) ;
		if ( i_out_of_date_ ) {
			*i_in_ = i;
			*i_out_ = index2multi_index2index_mod_shape( *dim_factor_out_, *shape_out_, *dim_factor_in_, *i_in_, *j_in_, *j_out_ ) ;
			i_out_of_date_ = false ;
		}
		else {
			index2multi_index2index_mod_shape_memory( *dim_factor_out_, *shape_out_, *dim_factor_in_, *shape_in_, *j_in_, *j_out_, i, *i_in_, *i_out_) ;
		}
	}

	void inc_lin_in ( ) const {
		if ( i_out_of_date_ ) { update_i () ; }
		auto i = *i_in_ + 1 ;
		if ( i == lin_size_ ) {
			++overflow_count_ ;
			i = 0 ;
		}
		index2multi_index2index_mod_shape_memory( *dim_factor_out_, *shape_out_, *dim_factor_in_, *shape_in_, *j_in_, *j_out_, i, *i_in_, *i_out_) ;
	}

	void dec_lin_in ( ) const {
		if ( i_out_of_date_ ) { update_i () ; }
		auto i = *i_in_ - 1 ;
		if ( i == -1 ) {
			--overflow_count_ ;
			i = lin_size_ - 1 ;
		}
		index2multi_index2index_mod_shape_memory( *dim_factor_out_, *shape_out_, *dim_factor_in_, *shape_in_, *j_in_, *j_out_, i, *i_in_, *i_out_) ;
	}

	void add_lin_in ( value_type const& j ) const {
		if ( i_out_of_date_ ) { update_i () ; }
		auto i = *i_in_ + j ;
		if ( i >= lin_size_ ) {
			std::ldiv_t dv{} ;
			dv = std::ldiv( long(i), long(lin_size_) ) ;
			i = dv.rem ;
			overflow_count_ += dv.quot ;
		}
		if ( i < 0 ) {
			std::ldiv_t dv{} ;
			dv = std::ldiv( - long(i), long(lin_size_) ) ;
			if ( dv.rem > 0 ) {
				i = - dv.rem + lin_size_ ;
			}
			overflow_count_ -= dv.quot + 1 ;
		}
		index2multi_index2index_mod_shape_memory( *dim_factor_out_, *shape_out_, *dim_factor_in_, *shape_in_, *j_in_, *j_out_, i, *i_in_, *i_out_) ;
	}

	value_type lin_in ( ) const {
		if ( i_out_of_date_ ) { update_i () ; }
		return *i_in_ ;
	}

	value_type lin_out ( ) const {
		if ( i_out_of_date_ ) { update_i () ; }
		return *i_out_ ;
	}

	index_type const& shape_in () const {
		return *shape_in_ ;
	}

	index_type const& shape_out () const {
		return *shape_out_ ;
	}

	template <typename J>
	void set_multi_in ( J const& j ) const {
		*j_in_ = j ;
		for ( size_type k = 0; k < j_in_->size(); ++k ) {
			assert( (*j_in_)[k] >= 0 && (*j_in_)[k] < (*shape_in_)[k] ) ;
			(*j_out_)[k] = (*j_in_)[k] % (*shape_out_)[k] ;
		}
		i_out_of_date_ = true ;
	}

	index_type const& multi_in ( ) const {
		return *j_in_ ;
	}

	index_type const& multi_out ( ) const {
		return *j_out_ ;
	}

	void reset( ) const {
		for ( ndims_type k = 0; k < shape_in_->size(); ++k ){
			(*j_in_)[k] = 0 ;
			(*j_out_)[k] = 0 ;
		}
		*i_in_ = 0 ;
		*i_out_ = 0 ;
		overflow_count_ = 0 ;
	}

public:
    template <typename I>
	auto operator[] ( I const& i ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( (*j_out_)[i] ) >::type {
		return (*j_out_)[i] ;
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
	auto operator() ( X const& x ) const -> typename std::enable_if< is<Vector, X>::value, decltype( (*j_out_)[x[0]] ) >::type {
		assert( x.size() == 1 ) ;
		return (*j_out_)[x[0]] ;
	}

	auto operator() ( std::initializer_list<size_type> const& x ) const -> decltype( (*j_out_)[*x.begin()] ) {
		assert( x.size() == 1 ) ;
		return (*j_out_)[*x.begin()] ;
	}

	template <typename = void>
	auto operator() ( std::vector<primitive_vector_wrapper<size_type>> const& s ) const -> decltype ( block_selection_of_vector( *this, *s.begin() ) ) {
		assert( s.size() == 1 ) ;
		return block_selection_of_vector( *this, *s.begin() ) ;
	}

} ;

template <typename V>
struct concept< tile_index<V>, typename std::enable_if< std::is_integral<V>::value >::type >
: Index
  {} ;

} // namespace glas3


#endif
