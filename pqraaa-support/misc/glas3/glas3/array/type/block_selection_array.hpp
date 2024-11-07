//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_block_selection_array_hpp
#define glas3_array_type_block_selection_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/primitive_vector_wrapper.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/primitive_dense_vector.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/assign.hpp>

#include <initializer_list>
#include <vector>
#include <cassert>
#include <type_traits>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

namespace glas3 {

template <typename T>
class block_selection_array {
	static_assert( is<Array, T>::value, "T should be an Array" ) ;

public:
	typedef T                               array_type ;
	typedef typename T::value_type          value_type ;
	typedef typename T::size_type           size_type ;
	typedef typename T::ndims_type          ndims_type ;
	typedef dense_vector<size_type>         shape_type ;

public:
	block_selection_array ( )
	: array_( boost::make_shared< array_type >( ) )
    , shape_( boost::make_shared< shape_type >( ) )
	, selection_arrays_( )
	, dim_factor_( boost::make_shared< shape_type >( ) )
	, size_( 0 )
	{}

	block_selection_array ( array_type const& array, std::vector<boost::shared_ptr<primitive_vector_wrapper<size_type>>> const& selection_arrays )
	: array_( boost::make_shared< array_type >( array.shallow_copy() ) )
	, shape_( boost::make_shared< shape_type >( no_init(), array.shape().size() ) )
	, selection_arrays_( boost::make_shared<std::vector<std::vector<boost::shared_ptr<primitive_vector_wrapper<size_type>>>>>( 1, selection_arrays ) )
	, dim_factor_( boost::make_shared< shape_type >( 1, array.shape().size() ) )
	, size_( 1 )
	{
		assert( array.shape().size() == selection_arrays.size() ) ;
	    for ( ndims_type k = 0; k < array.shape().size(); ++k ) {
			(*shape_)[k] = selection_arrays[k]->size() ;
			size_ *= (*shape_)[k] ;
			if ( k > 0 ) {
			    (*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ;
			}
		}
	}

	block_selection_array ( array_type const& array, std::vector<std::vector<boost::shared_ptr<primitive_vector_wrapper<size_type>>>> const& selection_arrays )
	: array_( boost::make_shared< array_type >( array.shallow_copy() ) )
	, shape_( boost::make_shared< shape_type >( no_init(), array.shape().size() ) )
	, selection_arrays_( boost::make_shared<std::vector<std::vector<boost::shared_ptr<primitive_vector_wrapper<size_type>>>>>( selection_arrays ) )
	, dim_factor_( boost::make_shared< shape_type >( 1, array.shape().size() ) )
	, size_( 1 )
	{
	    for ( ndims_type k = 0; k < array.shape().size(); ++k ) {
			(*shape_)[k] = selection_arrays.back()[k]->size() ;
			size_ *= (*shape_)[k] ;
			if ( k > 0 ) {
			    (*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ;
			}
		}
	}

public:
	// Copy constructor -> deep copy
	block_selection_array ( block_selection_array const& that )
	: array_( boost::make_shared< array_type >( *that.array_ ) )
	, shape_( boost::make_shared< shape_type >( *that.shape_in_ ) )
	, selection_arrays_( boost::make_shared<std::vector<std::vector<primitive_vector_wrapper<size_type>>>>( that.selection_arrays_-size() ) )
	, dim_factor_( boost::make_shared< shape_type >( *that.dim_factor_in_ ) )
	, size_( that.size_ )
	{
		for ( std::ptrdiff_t i = 0; i < selection_arrays_->size(); ++i ) {
			for ( auto selection_arrays_it = that.selection_arrays_->at(i)->begin(); selection_arrays_it != that.selection_arrays_->at(i)->end(); ++selection_arrays_it ) {
				selection_arrays_->at(i)->push_back( boost::make_shared<primitive_vector_wrapper<size_type>>( primitive_dense_vector<size_type>( *selection_arrays_it ) ) ) ;
			}
		}
	}

	// Move constructor
	block_selection_array ( block_selection_array && that ) = default ;

public:
	// Copy assignment -> deep copy
	block_selection_array& operator= ( block_selection_array const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	// Move assignment
	block_selection_array& operator= ( block_selection_array&& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	block_selection_array& operator= ( value_type const& value ) {
		assert( size() == 1 );
		(*this)[0] = value ;
		return *this;
	}

	block_selection_array& operator= ( std::initializer_list<value_type> const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, block_selection_array& >::type operator= ( E const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, block_selection_array& >::type operator= ( E const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

public:
	// Constructor used in shallow_copy
	block_selection_array ( boost::shared_ptr<array_type> array,
			boost::shared_ptr<shape_type> shape,
			boost::shared_ptr<std::vector<std::vector<boost::shared_ptr<primitive_vector_wrapper<size_type>>>>> selection_arrays,
			boost::shared_ptr<shape_type> dim_factor, size_type size )
	: array_( array )
	, shape_( shape )
	, selection_arrays_( selection_arrays )
	, dim_factor_( dim_factor )
	, size_( size )
	{}

private:
	boost::shared_ptr<array_type>                                                                          array_ ;
	boost::shared_ptr<shape_type>                                                                          shape_ ;
	boost::shared_ptr<std::vector<std::vector<boost::shared_ptr<primitive_vector_wrapper<size_type>>>>>    selection_arrays_ ;
	boost::shared_ptr<shape_type>                                                                          dim_factor_ ;
	size_type                                                                                              size_ ;

public:
	block_selection_array shallow_copy () const {
		return block_selection_array( array_, shape_, selection_arrays_, dim_factor_, size_ );
	}

public:
	shape_type shape () const {
		return *shape_ ;
	}

	size_type size () const {
		return size_ ;
	}

	size_type ndof () const {
		return array_->ndof() ;
	}

public:
	template <typename I>
	auto operator[] ( I const& i ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( (*array_)(shape_type()) ) >::type {
		assert( i >= 0 && i < size_ ) ;
		shape_type j( no_init(), shape_->size() ) ;
		index2multi_index( *dim_factor_, i, j ) ;
		for ( auto selection_arrays_it = selection_arrays_->rbegin(); selection_arrays_it != selection_arrays_->rend(); ++selection_arrays_it ) {
			for ( ndims_type k = 0; k < j.size(); ++k ) {
				j[k] = (*selection_arrays_it->at( k ))[j[k]] ;
			}
		}
		return (*array_)(j) ;
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
	auto operator() ( J const& j ) const -> typename std::enable_if< is<Vector, J>::value, decltype ( (*array_)(shape_type()) ) >::type {
		assert( shape_->size() == j.size() ) ;
		auto selection_arrays_it = selection_arrays_->rbegin() ;
		shape_type jj( no_init(), shape_->size() ) ;
		for ( ndims_type k = 0; k < j.size(); ++k ) {
			jj[k] = (*selection_arrays_it->at(k))[j[k]] ;
		}
		++selection_arrays_it ;
		for ( ; selection_arrays_it != selection_arrays_->rend(); ++selection_arrays_it ) {
			for ( ndims_type k = 0; k < j.size(); ++k ) {
				jj[k] = (*selection_arrays_it->at( k ))[jj[k]] ;
			}
		}
		return (*array_)(jj) ;
	}

	auto operator() ( std::initializer_list<size_type> const& j ) const -> decltype ( (*array_)(shape_type()) ) {
		assert( shape_->size() == j.size() ) ;
		auto selection_arrays_it = selection_arrays_->rbegin() ;
		shape_type jj( no_init(), shape_->size() ) ;
		ndims_type k = 0 ;
		for ( auto j_it = j.begin(); j_it != j.end(); ++j_it, ++k ) {
			jj[k] = (*selection_arrays_it->at(k))[*j_it] ;
		}
		++selection_arrays_it ;
		for ( ; selection_arrays_it != selection_arrays_->rend(); ++selection_arrays_it ) {
			for ( ndims_type k = 0; k < j.size(); ++k ) {
				jj[k] = (*selection_arrays_it->at( k ))[jj[k]] ;
			}
		}
		return (*array_)(jj) ;
	}

	template <typename = void>
	block_selection_array operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const {
		assert( shape_->size() == s.size() ) ;
		auto selection_arrays = *selection_arrays_ ;
		std::vector<boost::shared_ptr<primitive_vector_wrapper< size_type >>> dum( s.size() ) ;
		typename T::ndims_type k = 0 ;
		for ( auto s_it = s.begin(); s_it != s.end(); ++s_it, ++k ){
			dum[k] = boost::make_shared< primitive_vector_wrapper< size_type > >( s_it->shallow_copy() ) ;
		}
		selection_arrays.push_back( std::move( dum ) ) ;
		return block_selection_array ( *array_, selection_arrays ) ;
	}

} ;

template <typename T>
struct concept<block_selection_array<T>, typename std::enable_if< is<Array, T>::value && !is<DenseMatrix, T>::value
                                                  && !is<DenseVector, T>::value && !is<DenseScalar, T>::value >::type >
: DenseArray
  {};

template <typename T>
struct concept<block_selection_array<T>, typename std::enable_if< is<Matrix, T>::value >::type >
: DenseMatrix
  {};

template <typename T>
struct concept<block_selection_array<T>, typename std::enable_if< is<Vector, T>::value >::type >
: DenseVector
  {};

template <typename T>
struct concept<block_selection_array<T>, typename std::enable_if< is<Scalar, T>::value >::type >
: DenseScalar
  {};

} // namespace glas3


#endif
