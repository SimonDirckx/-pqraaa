//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_empty_array_hpp
#define glas3_array_type_empty_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <initializer_list>
#include <vector>
#include <type_traits>
#include <cassert>

namespace glas3 {

class empty_array {
public:
	typedef std::ptrdiff_t   value_type;
	typedef empty_array      shape_type ;
	typedef std::ptrdiff_t   size_type ;
	typedef std::ptrdiff_t   ndims_type ;

public:
	empty_array () {}

	// Copy constructor from std::initializer_list< value_type >
	empty_array ( std::initializer_list< value_type > data ) {
		assert( data.size() == 0 ) ;
	}

	// Copy constructor from Array
	template <typename E>
	empty_array ( E const& that ) {
		assert( that.size() == 0 ) ;
	}

public:
	// Copy constructor
	empty_array ( empty_array const& that ) = default ;

	// Move constructor
	empty_array ( empty_array && ) = default ;

public:
	// Copy assignment
	empty_array& operator= ( empty_array const& that ) {
		return *this ;
	}

	// Move assignment
	empty_array& operator= ( empty_array&& that ) = default;

	empty_array& operator= ( std::initializer_list< value_type > that ) {
		assert( that.size() == 0 );
		return *this ;
	}

	template <typename E>
	empty_array& operator= ( E const& that ) {
		assert( that.size() == 0 ) ;
		return *this ;
	}

public:
	empty_array shallow_copy () const {
		return empty_array () ;
	}

public:
	empty_array shape () const {
		return empty_array() ;
	}

	size_type size () const {
		return 0 ;
	}

	size_type ndof () const {
		return 0 ;
	}

public:
	template <typename I>
	typename std::enable_if< std::is_integral<I>::value, value_type >::type
	operator[] ( I const& i ) const {
		assert( false ) ;
		return 0 ;
	}

	template <typename X>
	typename std::enable_if< is<Array, X>::value, empty_array >::type
	operator[] ( X const& x ) const {
		assert( x.size() == 0 ) ;
		return this->shallow_copy() ;
	}

	template <typename I>
	typename std::enable_if< std::is_integral<I>::value, empty_array >::type
	operator[] ( std::initializer_list<I> x ) const {
		assert( x.size() == 0 ) ;
		return this->shallow_copy() ;
	}

	template <typename X>
	typename std::enable_if< is<Vector, X>::value, value_type >::type
	operator() ( X const& x ) const {
		assert( false ) ;
		return 0 ;
	}

	empty_array
	operator() ( std::initializer_list<size_type> j ) const {
		assert( j.size() == 0 ) ;
		return this->shallow_copy() ;
	}

} ;

template <>
struct concept< empty_array >
: DenseArray
  {};

} // namespace glas3


#endif
