#include <glas3/array/type/empty_array.hpp>
#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>

#include <glas3/array/dense_array/algorithm/ops.hpp>
#include <glas3/array/dense_array/algorithm/permute.hpp>
#include <glas3/array/dense_array/algorithm/ttt.hpp>
#include <glas3/array/algorithm/ndims.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

 	std::cout << "ndims operation" << '\n' ;
 	glas3::empty_array e ;
 	glas3::dense_scalar<int> s = 1 ;
 	glas3::dense_vector<int> v = {1, 2, 3, 4, 5, 6} ;
 	glas3::dense_matrix<int> m({1, 2, 3, 4, 5, 6}, {2, 3}) ;
 	glas3::dense_array<int> a({1, 2, 3, 4, 5, 6}, {2, 3, 1}) ;
    auto b = ( a * ( v + a ) - v ) / v ;
    auto p = glas3::permute( a, {2, 0, 1} ) ;
    auto t = glas3::ttt( v, a ) ;

    if ( glas3::ndims(e) != e.shape().size() || glas3::ndims(v) != v.shape().size() || glas3::ndims(m) != m.shape().size()
    		|| glas3::ndims(s) != s.shape().size() || glas3::ndims(a) != a.shape().size()
    		|| glas3::ndims(b) != b.shape().size() || glas3::ndims(p) != p.shape().size()
    		|| glas3::ndims(t) != t.shape().size() ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
