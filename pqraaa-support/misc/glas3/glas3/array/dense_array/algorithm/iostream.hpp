//  (C) Copyright sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_iostream_hpp
#define glas3_array_dense_array_algorithm_iostream_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>
#include <glas3/array/dense_array/concept/contiguous_dense_array.hpp>

#include <iostream>
#include <type_traits>

namespace glas3 {

template <typename V>
typename std::enable_if< is<DenseArray, V>::value, std::ostream& >::type operator<< ( std::ostream& s, V const& v ) {

	s << "concept: " ;

	if ( is<DenseArray, V>::value ) {
		if ( is<ContiguousDenseArray, V>::value ) {
			if ( is<ContiguousDenseScalar, V>::value ) { s << "ContiguousDenseScalar" << '\n' ; }
			else if ( is<ContiguousDenseVector, V>::value ) { s << "ContiguousDenseVector" << '\n' ; }
			else if ( is<ContiguousDenseMatrix, V>::value ) { s << "ContiguousDenseMatrix" << '\n' ; }
			else { s << "ContiguousDenseArray" << '\n' ; }
		}
		else {
			if ( is<DenseScalar, V>::value ) { s << "DenseScalar" << '\n' ; }
			else if ( is<DenseVector, V>::value ) { s << "DenseVector" << '\n' ; }
			else if ( is<DenseMatrix, V>::value ) { s << "DenseMatrix" << '\n' ; }
			else { s << "DenseArray" << '\n' ; }
		}
	}

	s << "shape: { " ;
	if ( v.shape().size() > 0 ) s << v.shape()[0] ;
	for ( typename V::ndims_type i = 1; i < v.shape().size(); ++i) { s << ", " << v.shape()[i] ; }
	s << " }\n" ;

	s << "data: { " ;
	if ( v.size() > 0 ) s << v[0] ;
	for ( typename V::size_type i = 1; i < v.size(); ++i) { s << ", " << v[i] ; }
	s << " }\n" ;

	return s ;
}

} // namespace glas3

#endif
