//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_block_selection_vector_hpp
#define glas3_array_type_block_selection_vector_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/concept/container.hpp>
#include <glas3/concept/random_access_container.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/primitive_dense_vector.hpp>
#include <glas3/array/type/primitive_vector_wrapper.hpp>

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
class block_selection_vector {
	static_assert( is<Vector, T>::value, "T should be a Vector" ) ;

public:
	typedef T                               vector_type ;
	typedef typename T::value_type          value_type ;
	typedef typename T::size_type           size_type ;
	typedef typename T::ndims_type          ndims_type ;
	typedef dense_scalar<size_type>         shape_type ;

public:
	block_selection_vector ( )
	: vector_( boost::make_shared< vector_type >( ) )
    , selection_vectors_( )
    , shape_( boost::make_shared< shape_type >( ) )
	, size_( 0 )
	{}

	block_selection_vector ( vector_type const& vector, primitive_vector_wrapper<size_type> const& selection_vector )
	: vector_( boost::make_shared< vector_type >( vector.shallow_copy() ) )
    , selection_vectors_( boost::make_shared<std::vector<primitive_vector_wrapper<size_type>>>( 1, selection_vector ) )
	, shape_( boost::make_shared< shape_type >( selection_vector.size() ) )
	, size_( selection_vector.size() )
	{}

	block_selection_vector ( vector_type const& vector, std::vector<primitive_vector_wrapper<size_type>> const& selection_vectors )
	: vector_( boost::make_shared< vector_type >( vector.shallow_copy() ) )
	, selection_vectors_( boost::make_shared<std::vector<primitive_vector_wrapper<size_type>>>( selection_vectors ) )
	, shape_( boost::make_shared< shape_type >( selection_vectors.back().size() ) )
	, size_( selection_vectors.back().size() )
	{}

public:
	// Copy constructor -> deep copy
	block_selection_vector ( block_selection_vector const& that )
	: vector_( boost::make_shared< vector_type >( *that.vector_ ) )
	, selection_vectors_( boost::make_shared<std::vector<primitive_vector_wrapper<size_type>>>( ) )
	, shape_( boost::make_shared< shape_type >( *that.shape_ ) )
	, size_( that.size_ )
	{
		for ( auto selection_vectors_it = that.selection_vectors_->begin(); selection_vectors_it != that.selection_vectors_->end(); ++selection_vectors_it ) {
			selection_vectors_->push_back( primitive_vector_wrapper<size_type>( primitive_dense_vector<size_type>( *selection_vectors_it ) ) ) ;
		}
	}

	// Move constructor
	block_selection_vector ( block_selection_vector && that ) = default ;

public:
	// Copy assignment -> deep copy
	block_selection_vector& operator= ( block_selection_vector const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	// Move assignment
	block_selection_vector& operator= ( block_selection_vector&& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

	block_selection_vector& operator= ( value_type const& value ) {
		assert( size() == 1 );
		(*this)[0] = value ;
		return *this;
	}

	block_selection_vector& operator= ( std::initializer_list<value_type> const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Container, E>::value, block_selection_vector& >::type operator= ( E const& that ) {
		assert( that.size() == size() );
		size_type i = 0 ;
		for ( auto that_it = that.begin(); that_it != that.end(); ++that_it, ++i ) {
			(*this)[i] = *that_it ;
		}
		return *this;
	}

	template <typename E>
	typename std::enable_if< is<Array, E>::value, block_selection_vector& >::type operator= ( E const& that ) {
		assert( that.size() == size() ) ;
		assign( *this, that ) ;
		return *this ;
	}

public:
	// Constructor used in shallow_copy
	block_selection_vector ( boost::shared_ptr<vector_type> vector,
			std::vector<primitive_vector_wrapper<size_type>> selection_vectors,
			boost::shared_ptr<shape_type> shape, size_type size )
	: vector_( vector )
    , selection_vectors_( selection_vectors )
	, shape_( shape )
	, size_( size )
	{}

private:
	boost::shared_ptr<vector_type>                                         vector_ ;
	boost::shared_ptr<std::vector<primitive_vector_wrapper<size_type>>>    selection_vectors_ ;
	boost::shared_ptr<shape_type>                                          shape_ ;
	size_type                                                              size_ ;

public:
	block_selection_vector shallow_copy () const {
		return block_selection_vector( vector_, selection_vectors_, shape_, size_ );
	}

public:
	shape_type shape () const {
		return *shape_ ;
	}

	size_type size () const {
		return size_ ;
	}

	size_type ndof () const {
		size_type n = vector_->ndof() ;
		for ( auto selection_vectors_it = selection_vectors_->begin(); selection_vectors_it != selection_vectors_->end(); ++selection_vectors_it ) {
			n += selection_vectors_it->ndof() ;
		}
		return n ;
	}

public:
	template <typename I>
	auto operator[] ( I const& i ) const -> typename std::enable_if< std::is_integral<I>::value, decltype( (*vector_)[size_type()] ) >::type {
		assert( i >= 0 && i < size_ ) ;
		auto selection_vectors_it = selection_vectors_->rbegin() ;
		size_type ii = (*selection_vectors_it)[i] ; ++selection_vectors_it ;
		for ( ; selection_vectors_it != selection_vectors_->rend(); ++selection_vectors_it ) {
			ii = (*selection_vectors_it)[ii] ;
		}
		return (*vector_)[ii] ;
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
	auto operator() ( J const& j ) const -> typename std::enable_if< is<Vector, J>::value, decltype ( (*vector_)[size_type()] ) >::type {
		assert( j.size() == 1 ) ;
		auto selection_vectors_it = selection_vectors_->rbegin() ;
		size_type ii = (*selection_vectors_it)[j[0]] ; ++selection_vectors_it ;
		for ( ; selection_vectors_it != selection_vectors_->rend(); ++selection_vectors_it ) {
			ii = (*selection_vectors_it)[ii] ;
		}
		return (*vector_)[ii] ;
	}

	auto operator() ( std::initializer_list<size_type> const& j ) const -> decltype ( (*vector_)[size_type()] ) {
		assert( j.size() == 1 ) ;
		auto selection_vectors_it = selection_vectors_->rbegin() ;
		size_type ii = (*selection_vectors_it)[*j.begin()] ; ++selection_vectors_it ;
		for ( ; selection_vectors_it != selection_vectors_->rend(); ++selection_vectors_it ) {
			ii = (*selection_vectors_it)[ii] ;
		}
		return (*vector_)[ii] ;
	}

	template <typename = void>
	block_selection_vector operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const {
		assert( s.size() == 1 ) ;
		auto selection_vectors = *selection_vectors_ ;
		selection_vectors.push_back( *s.begin() ) ;
		return block_selection_vector ( *vector_, selection_vectors ) ;
	}

} ;

template <typename T>
struct concept<block_selection_vector<T>, typename std::enable_if< is<Vector, T>::value >::type >
: DenseVector
  {};

} // namespace glas3


#endif
