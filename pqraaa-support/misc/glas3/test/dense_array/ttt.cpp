#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ttt.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k ;

	std::cout << "ttt operation (shallow copy): outer product" << '\n' ;
	std::cout << "auto ttt_va = glas3::ttt( glas3::dense_vector<int>({2, 3}), glas3::dense_array<int>({5, 7, 11, 13}, {2, 2}) ) ;" << '\n' ;
 	glas3::dense_vector<int> v1 = {100, 3} ;
 	glas3::dense_array<int> a1({5, 7, 11, 13}, {2, 2}) ;
    auto ttt_v1a1_2 = glas3::ttt( v1, a1 ) ;
    if ( ttt_v1a1_2.size() != 8 || ttt_v1a1_2.ndof() != 2 + 4 || ttt_v1a1_2.shape().size() != 3 || ttt_v1a1_2.shape()[0] != 2 || ttt_v1a1_2.shape()[1] != 2 || ttt_v1a1_2.shape()[2] != 2 || ttt_v1a1_2.shape().shape().size() != 1 ) return 1 ;
    v1[0] -= 1 ;
    for ( i = 0; i < 2; ++i ) {
    	for ( j = 0; j < 2; ++j ) {
    		for ( k = 0; k < 2; ++k ) {
       			if ( ttt_v1a1_2({i, j, k}) != v1({i}) * a1({j, k}) || ttt_v1a1_2[i + 2 * j + 2 * 2 * k] != v1[i] * a1[j + 2 * k] ) return 1 ;
    		}
    	}
    }

    std::cout << "ttt operation (shallow copy): outer product" << '\n' ;
    std::cout << "auto ttt_va = glas3::ttt( glas3::dense_scalar<int>(2), glas3::dense_array<int>({5, 7, 11, 13}, {2, 2}) ) ;" << '\n' ;
 	glas3::dense_scalar<int> s1 = 2 ;
 	glas3::dense_array<int> a2({5, 7, 11, 13}, {2, 2}) ;
    auto ttt_s1a2 = glas3::ttt( s1, a2 ) ;
    if ( ttt_s1a2.size() != 4 || ttt_s1a2.ndof() != 1 + 4 || ttt_s1a2.shape().size() != 2 || ttt_s1a2.shape()[0] != 2 || ttt_s1a2.shape()[1] != 2 || ttt_s1a2.shape().shape().size() != 1 ) return 1 ;
    s1[0] -= 1 ;
    for ( j = 0; j < 2; ++j ) {
    	for ( k = 0; k < 2; ++k ) {
    		if ( ttt_s1a2({j, k}) != s1[0] * a2({j, k}) || ttt_s1a2[j + 2 * k] != s1[0] * a2[j + 2 * k] ) return 1 ;
    	}
    }

	std::cout << "ttt operation (shallow copy): contraction along some dimensions" << '\n' ;
	std::cout << "auto ttt_va = glas3::ttt( glas3::dense_vector<int>({2, 3}), glas3::dense_array<int>({5, 7, 11, 13}, {2, 2}), {0}, {1} ) ;" << '\n' ;
 	glas3::dense_vector<int> v2 = {2, 3} ;
 	glas3::dense_array<int> a3({5, 7, 11, 13}, {2, 2}) ;
    auto ttt_v2a3 = glas3::ttt( v2, a3, {0}, {1} ) ;
    if ( ttt_v2a3.size() != 2 || ttt_v2a3.ndof() != 2 + 4 || ttt_v2a3.shape().size() != 1 || ttt_v2a3.shape()[0] != 2 || ttt_v2a3.shape().shape().size() != 1 ) return 1 ;
    v2[0] -= 1 ;
    for ( j = 0; j < 2; ++j ) {
    	int dum1 = 0;
    	for ( k = 0; k < 2; ++k ) {
        	dum1 += v2({k}) * a3({j, k}) ;
   		}
   		if ( ttt_v2a3({j}) != dum1 || ttt_v2a3[j] != dum1 ) return 1 ;
    }

	std::cout << "ttt operation (shallow copy): inner product" << '\n' ;
	std::cout << "auto ttt_va = glas3::ttt( glas3::dense_array<int>({2, 3, 5, 7}, {2, 2}), glas3::dense_array<int>({11, 13, 17, 19}, {2, 2}), {0, 1}, {1, 0} ) ;" << '\n' ;
 	glas3::dense_array<int> a4({2, 3, 5, 7}, {2, 2}) ;
 	glas3::dense_array<int> a5({11, 13, 17, 19}, {2, 2}) ;
    auto ttt_a4a5 = glas3::ttt( a4, a5, {0, 1}, {1, 0} ) ;
    if ( ttt_a4a5.size() != 1 || ttt_a4a5.ndof() != 4 + 4 || ttt_a4a5.shape().size() != 0 || ttt_a4a5.shape().shape().size() != 1 ) return 1 ;
    a4[0] -= 1 ;
    int dum1 = 0;
    for ( j = 0; j < 2; ++j ) {
    	for ( k = 0; k < 2; ++k ) {
        	dum1 += a4({j, k}) * a5({k, j}) ;
   		}
    }
    if ( ttt_a4a5[0] != dum1 ) return 1 ;

	std::cout << '\n' ;

	return 0;

}
