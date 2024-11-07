//	  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_entrywise_operation_hpp
#define glas3_array_dense_array_type_entrywise_operation_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/type/entrywise_operation.hpp>
#include <glas3/array/type/primitive_vector_wrapper.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>
#include <glas3/array/dense_array/algorithm/block_selection.hpp>
#include <glas3/array/algorithm/indexing.hpp>

#include <initializer_list>
#include <vector>
#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>
#include <functional>

namespace glas3 {

template < typename Op, typename X >
class unary_operation< Op, X, typename std::enable_if< is< DenseArray, X >::value >::type > {
public:
	typedef typename X::shape_type          shape_type ;
	typedef typename X::size_type           size_type ;
	typedef typename X::ndims_type          ndims_type ;
	typedef typename Op::result_type        value_type ;

public:
	unary_operation ( )
    : op_ ( )
	, x_( boost::make_shared< X >( ) )
    {}

	unary_operation ( Op op, X const& x )
    : op_ ( op )
	, x_( boost::make_shared< X >( x.shallow_copy() ) )
    {}

public:
    // Copy constructor -> deep copy
	unary_operation ( unary_operation const& that )
	: op_ ( that.op_ )
	, x_( boost::make_shared< X >( *that.x_ ) )
	{}

    // Move constructor
	unary_operation ( unary_operation && that ) = default;

public:
    // Copy assignment -> do not allow
	unary_operation& operator= ( unary_operation const& that ) = delete ; // not assignable

    // Move assignment
	unary_operation& operator= ( unary_operation&& that ) = delete ;

private:
    // Constructor used in shallow_copy
	unary_operation ( Op op, boost::shared_ptr<X> x )
    : op_ ( op )
	, x_( x )
    {}

private:
	Op                       op_ ;
    boost::shared_ptr< X >   x_;

public:
	unary_operation shallow_copy () const {
       	return unary_operation( op_, x_ );
    }

public:
	shape_type const& shape () const {
		return x_->shape() ;
	}

	size_type size () const {
		return x_->size() ;
	}

	size_type ndof () const {
		return x_->ndof() ;
	}

public:
    template <typename I>
	typename std::enable_if< std::is_integral<I>::value, value_type >::type operator[] ( I const& i ) const {
		return op_( x_->operator[](i) ) ;
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
    typename std::enable_if< is<Vector, S>::value, value_type >::type operator() ( S const& x ) const {
    	assert( x.size() == x_->shape().size() ) ;
    	return op_( x_->operator[]( multi_index2index_initializer_list( x_->shape(), x ) ) ) ;
    }

	value_type operator() ( std::initializer_list<size_type> const& x ) const {
		assert( x.size() == x_->shape().size() ) ;
		return op_( x_->operator[]( multi_index2index_initializer_list( x_->shape(), primitive_vector_wrapper<size_type>( x ) ) ) ) ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

} ;

template <typename Op, typename X>
struct concept< unary_operation<Op, X>, typename std::enable_if< is<DenseArray, X>::value
                                                            && ! is<DenseMatrix, X>::value
                                                            && ! is<DenseVector, X>::value
                                                            && ! is<DenseScalar, X>::value >::type >
: DenseArray
  {} ;

template <typename Op, typename X>
struct concept< unary_operation<Op, X>, typename std::enable_if< is<DenseMatrix, X>::value >::type >
: DenseMatrix
  {} ;

template <typename Op, typename X>
struct concept< unary_operation<Op, X>, typename std::enable_if< is<DenseVector, X>::value >::type >
: DenseVector
  {} ;

template <typename Op, typename X>
struct concept< unary_operation<Op, X>, typename std::enable_if< is<DenseScalar, X>::value >::type >
: DenseScalar
  {} ;

template < typename Op, typename X, typename Y >
class binary_operation< Op, X, Y, typename std::enable_if< is<DenseArray, X>::value && is<DenseArray, Y>::value >::type > {
public:
	typedef typename X::shape_type        shape_type ;
	typedef typename X::size_type         size_type ;
	typedef typename X::ndims_type        ndims_type ;
	typedef typename Op::result_type      value_type ;

public:
	binary_operation ( Op op, X const& x, Y const& y )
    : op_ ( op )
	, x_( boost::make_shared< X >( x.shallow_copy() ) )
    , y_( boost::make_shared< Y >( y.shallow_copy() ) )
    {
		assert( x.size() == y.size() );
    }

public:
    // Copy constructor -> deep copy
	binary_operation ( binary_operation const& that )
    : op_ ( that.op_ )
    , x_( boost::make_shared< X >( *that.x_ ) )
    , y_( boost::make_shared< Y >( *that.y_ ) )
    {}

    // Move constructor
	binary_operation ( binary_operation && that ) = default;

public:
    // Copy assignment -> do not allow
	binary_operation& operator= ( binary_operation const& that ) = delete ; // not assignable

    // Move assignment
	binary_operation& operator= ( binary_operation&& that ) = delete ;

private:
    // Constructor used in shallow_copy
	explicit binary_operation ( Op op, boost::shared_ptr<X> x, boost::shared_ptr<Y> y )
    : op_ ( op )
	, x_( x )
    , y_( y )
    {}

private:
	Op                       op_ ;
	boost::shared_ptr< X >   x_;
	boost::shared_ptr< Y >   y_;

public:
	binary_operation shallow_copy () const {
       	return binary_operation( op_, x_, y_ );
    }

public:
	shape_type const& shape () const {
		return x_->shape() ;
	}

	size_type size () const {
		return x_->size() ;
	}

	size_type ndof () const {
		return x_->ndof() + y_->ndof() ;
	}

public:
    template <typename I>
    typename std::enable_if< std::is_integral<I>::value, value_type >::type operator[] ( I const& i ) const {
        //return op_( (*x_)[i], (*y_)[i] ) ;
    	return op_( x_->operator[](i), y_->operator[](i) ) ;
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
    typename std::enable_if< is<Vector, S>::value, value_type >::type operator() ( S const& x ) const {
    	assert( x.size() == x_->shape().size() ) ;
    	auto i = multi_index2index_initializer_list( x_->shape(), x ) ;
    	return op_( x_->operator[](i), y_->operator[](i) ) ;
    }

	value_type operator() ( std::initializer_list<size_type> const& x ) const {
    	assert( x.size() == x_->shape().size() ) ;
    	auto i = multi_index2index_initializer_list( x_->shape(), primitive_vector_wrapper<size_type>( x ) ) ;
    	return op_( x_->operator[](i), y_->operator[](i) ) ;
	}

	template <typename = void>
	auto operator() ( std::vector< primitive_vector_wrapper< size_type > > const& s ) const -> decltype ( block_selection( *this, s ) ) {
		return block_selection( *this, s ) ;
	}

} ;

template < typename Op, typename X, typename Y >
struct concept< binary_operation< Op, X, Y >, typename std::enable_if< is<DenseArray, X>::value
                                                            && ! is<DenseMatrix, X>::value
                                                            && ! is<DenseVector, X>::value
                                                            && ! is<DenseScalar, X>::value && is< DenseArray, Y >::value >::type >
: DenseArray
  {} ;

template < typename Op, typename X, typename Y >
struct concept< binary_operation< Op, X, Y >, typename std::enable_if< is<DenseMatrix, X>::value && is< DenseArray, Y >::value >::type >
: DenseMatrix
  {} ;

template < typename Op, typename X, typename Y >
struct concept< binary_operation< Op, X, Y >, typename std::enable_if< is<DenseVector, X>::value && is< DenseArray, Y >::value >::type >
: DenseVector
  {} ;

template < typename Op, typename X, typename Y >
struct concept< binary_operation< Op, X, Y >, typename std::enable_if< is<DenseScalar, X>::value && is< DenseArray, Y >::value >::type >
: DenseScalar
  {} ;

//template <typename Op, typename X, typename... Y>
//class entrywise_operation
//{
//	static_assert(is<DenseArray, X>::value, "all operands should be of DenseArray type"); //&& is<DenseArray, Y...>::value
//
//private:
//	template<int ...> struct seq {};
//	template<int N, int ...S> struct gens : gens<N-1, N-1, S...> {};
//	template<int ...S> struct gens<0, S...>{ typedef seq<S...> type; };
//
//public:
//	entrywise_operation( X const& x, Y const&... y )
//    : x_( &x )
//    , y_( &y... )
//    {}
//
//public:
//	typedef typename X::shape_type                                                     shape_type ;
//	typedef typename X::size_type                                                      size_type ;
//	typedef typename X::ndims_type                                                     ndims_type ;
//	typedef decltype( Op() ( typename X::value_type(), typename Y::value_type()... ) ) value_type ;
//
//	shape_type const& shape() const {
//		return x_->shape() ;
//	}
//	size_type size() const {
//		return x_->size() ;
//	}
//
//	value_type operator[] ( size_type i) const { return callFunc(i, typename gens<sizeof...(Y)>::type()) ; }
//	//auto operator() ( size_type i) const -> decltype( this->op_( this->scalar_, this->vector_(0) ) )  { return op_( scalar_, vector_(i) ) ; }
//
//private:
//
//	template<int ...S>
//	value_type callFunc( size_type i, seq<S...> ) { return Op()( (*x_)[i], (*std::get<S>( y_ ))[i] ... ); }
//
//private:
//	X const* x_;
//	std::tuple<Y const*...> y_;
//} ;
//
//template <typename Op, typename X, typename... Y>
//struct concept< entrywise_operation<Op, X, Y...>, typename std::enable_if< is<DenseArray, X>::value>::type > //&& is<DenseArray, Y...>::value
//: DenseArray
//  {} ;

} // namespace glas3

#endif
