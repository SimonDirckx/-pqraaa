#include <glas3/array/array_network/type/CP_tensor.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>

#include <glas3/array/array_network/algorithm/full.hpp>
#include <glas3/array/array_network/algorithm/fill.hpp>

#include <glas3/dense_array.hpp>

#include <iostream>

int main() {

 	std::cout << '\n';

    std::cout << "tensor_train construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::tensor_train<glas3::dense_array<double>> t(0, {5, 4, 3}, {4, 3}) ;" << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s( {5, 4, 3} ) ;
    glas3::dense_vector<std::ptrdiff_t> r( {4, 3} ) ;
    glas3::tensor_train<glas3::dense_array<double>> t(0, s, r) ;
    glas3::fill( t, 5 ) ;

    for ( auto nodes2arrays_it: t.nodes2arrays() ) {
    	auto node = *nodes2arrays_it.second ;
    	if ( glas3::norm_inf( glas3::full( node ) - glas3::constant_array<double>( 5, node.shape() ) ) != 0 ) return 1 ;
    }

    std::cout << '\n' ;

	return 0;

}
