#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ttt_1D.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::ptrdiff_t i, j, k, l, m ;

	std::cout << "ttt_1D operation (shallow copy)" << '\n' ;
	std::cout << "auto ttt_va = glas3::ttt_1D( glas3::dense_vector<int>({2, 3}), glas3::dense_array<int>({5, 7, 11, 13}, {2, 2}), 0 ) ;" << '\n' ;
 	glas3::dense_vector<int> v1 = {2, 3} ;
 	glas3::dense_array<int> a1({5, 7, 11, 13}, {2, 2}) ;
    auto ttt_v1a1_1 = glas3::ttt_1D( v1, a1, 0, 0 ) ;
    auto ttt_v1a1_2 = glas3::ttt_1D( a1, v1, 1, 0 ) ;
    auto ttt_v1a1_3 = glas3::ttt_1D( a1, v1 ) ;
    if ( ttt_v1a1_1.size() != 2 || ttt_v1a1_1.ndof() != v1.ndof() + a1.ndof() || ttt_v1a1_1.shape().size() != 1 || ttt_v1a1_1.shape()[0] != 2 || ttt_v1a1_1.shape().shape().size() != 1 ) return 1 ;
    if ( ttt_v1a1_2.size() != 2 || ttt_v1a1_2.ndof() != v1.ndof() + a1.ndof() || ttt_v1a1_2.shape().size() != 1 || ttt_v1a1_2.shape()[0] != 2 || ttt_v1a1_2.shape().shape().size() != 1 ) return 1 ;
    if ( ttt_v1a1_3.size() != 2 || ttt_v1a1_3.ndof() != v1.ndof() + a1.ndof() || ttt_v1a1_3.shape().size() != 1 || ttt_v1a1_3.shape()[0] != 2 || ttt_v1a1_3.shape().shape().size() != 1 ) return 1 ;
    v1[0] -= 1 ;
    for ( j = 0; j < 2; ++j ) {
    	int dum1 = 0;
    	int dum2 = 0;
    	for ( k = 0; k < 2; ++k ) {
        	dum1 += v1({k}) * a1({k, j}) ;
        	dum2 += v1({k}) * a1({j, k}) ;
   		}
   		if ( ttt_v1a1_1({j}) != dum1 || ttt_v1a1_1[j] != dum1 ) return 1 ;
   		if ( ttt_v1a1_2({j}) != dum2 || ttt_v1a1_2[j] != dum2 ) return 1 ;
   		if ( ttt_v1a1_3({j}) != dum2 || ttt_v1a1_3[j] != dum2 ) return 1 ;
    }

	std::cout << "ttt operation (shallow copy): inner product" << '\n' ;
	std::cout << "auto ttt_va = glas3::ttt( glas3::dense_array<int>({2, 3, 5, 7}, {2, 2}), glas3::dense_array<int>({11, 13, 17, 19}, {2, 2}), {0, 1}, {1, 0} ) ;" << '\n' ;
 	glas3::dense_array<int> a4({2, 3, 5, 7, 11, 13}, {1, 2, 3}) ;
 	glas3::dense_array<int> a5({17, 19, 23, 29, 31, 37}, {3, 2, 1}) ;
    auto ttt_a4a5 = glas3::ttt_1D( a4, a5, 1, 1 ) ;
    if ( ttt_a4a5.size() != 9 || ttt_a4a5.ndof() != a4.ndof() + a5.ndof() || ttt_a4a5.shape().size() != 4
    		 || ttt_a4a5.shape()[0] != 1 || ttt_a4a5.shape()[1] != 3 || ttt_a4a5.shape()[2] != 1 || ttt_a4a5.shape()[3] != 3
    		 || ttt_a4a5.shape().shape().size() != 1 ) return 1 ;
    a4[0] -= 1 ;
    for ( i = 0; i < 1; ++i ) {
    	for ( j = 0; j < 3; ++j ) {
    		for ( k = 0; k < 1; ++k ) {
    			for ( l = 0; l < 3; ++l ) {
    				int dum1 = 0 ;
    				for ( m = 0; m < 2; ++m ) {
    					dum1 += a4({i, m, l}) * a5({j, m, k}) ;
    				}
    				if ( ttt_a4a5({i, j, k, l}) != dum1 || ttt_a4a5[i + 1 * j + 1 * 3 * k + 1 * 3 * 1 * l] != dum1 ) return 1 ;
    			}
    		}
   		}
    }

	std::cout << '\n' ;

	return 0;

}
