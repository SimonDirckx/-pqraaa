#ifndef glas3_array_dense_array_algorithm_ops_assign_hpp
#define glas3_array_dense_array_algorithm_ops_assign_hpp

#include <glas3/concept/is.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <type_traits>
#include <cassert>

namespace glas3 {

template <typename V, typename E>
typename std::enable_if< is<DenseArray, V>::value && is<DenseArray, E>::value, V const& >::type operator+=( V const& v, E const& e ) {
	for (typename V::size_type i = 0; i < v.size(); ++i) {
		v[i] += e[i] ;
	}
	return v ;
}

template <typename V, typename E>
typename std::enable_if< is<DenseArray,V>::value && is<DenseArray,E>::value, V const& >::type operator-=( V const& v, E const& e ) {
	assert( v.size() == e.size() ) ;
	for (typename V::size_type i = 0; i < v.size(); ++i) {
		v[i] -= e[i] ;
	}
	return v ;
}

template <typename V, typename E>
typename std::enable_if< is<DenseArray,V>::value && is<DenseArray,E>::value, V const& >::type operator*=( V const& v, E const& e ) {
	assert( v.size() == e.size() ) ;
	for (typename V::size_type i = 0; i < v.size(); ++i) {
		v[i] *= e[i] ;
	}
	return v ;
}

template <typename V, typename E>
typename std::enable_if< is<DenseArray,V>::value && is<DenseArray,E>::value, V const& >::type operator/=( V const& v, E const& e ) {
	assert( v.size() == e.size() ) ;
	for (typename V::size_type i = 0; i < v.size(); ++i) {
		v[i] /= e[i] ;
	}
	return v ;
}

template <typename V, typename E>
typename std::enable_if< is<DenseArray,V>::value && is<DenseArray,E>::value, V const& >::type operator%( V const& v, E const& e ) {
	assert( v.size() == e.size() ) ;
	for (typename V::size_type i = 0; i < v.size(); ++i) {
		v[i] %= e[i] ;
	}
	return v ;
}

} // namespace glas3

#endif
