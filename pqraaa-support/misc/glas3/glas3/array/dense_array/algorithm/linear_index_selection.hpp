//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_linear_index_selection_hpp
#define glas3_array_dense_array_algorithm_linear_index_selection_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/type/indirect_array.hpp>
#include <glas3/array/dense_array/container/primitive_dense_vector.hpp>

#include <glas3/array/algorithm/linear_index_selection.hpp>

#include <vector>
#include <initializer_list>
#include <type_traits>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

namespace glas3 {

template < typename T, typename S >
typename std::enable_if< is< DenseArray, T >::value && is< DenseArray, S >::value, indirect_array< T, S > >::type
linear_index_selection( T const& t, S const& s ) {
    return indirect_array< T, S >( t, s ) ;
}

template < typename T, typename V >
typename std::enable_if< is< DenseArray, T >::value && std::is_integral<V>::value, indirect_array< T, primitive_dense_vector<V> > >::type
linear_index_selection( T const& t, std::initializer_list<V> const& s ) {
    return indirect_array< T, primitive_dense_vector<V> >( t, primitive_dense_vector<V> ( s ) ) ;
}

} // namespace glas3

#endif
