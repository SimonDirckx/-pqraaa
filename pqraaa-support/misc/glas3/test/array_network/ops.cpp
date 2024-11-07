#include <glas3/array/array_network/type/CP_tensor.hpp>
#include <glas3/array/array_network/type/Tucker_tensor.hpp>
#include <glas3/array/array_network/type/tensor_train.hpp>

#include <glas3/array/array_network/algorithm/randomize.hpp>
#include <glas3/array/array_network/algorithm/full.hpp>
#include <glas3/array/array_network/algorithm/ttt.hpp>
#include <glas3/array/array_network/algorithm/ops.hpp>
#include <glas3/array/array_network/algorithm/iostream.hpp>

#include <glas3/dense_array.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <iostream>

int main() {

 	std::cout << '\n';

 	typedef typename boost::uniform_real<double>                        distribution_type ;
 	typedef typename boost::lagged_fibonacci607                         engine_type ;
 	typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
 	std::random_device rd ;
 	generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

    std::cout << "CP tensor construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::CP_tensor<glas3::dense_array<double>> t1( 0, {3, 4, 5}, 3 ) ;" << '\n' << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s1( {3, 4, 5} ) ;
    std::ptrdiff_t r1( 3 ) ;
    glas3::CP_tensor<glas3::dense_array<double>> t1( 0, s1, r1 ) ;
    glas3::randomize( t1, g ) ;


    std::cout << "CP tensor construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::CP_tensor<glas3::dense_array<double>> t11( 0, {3, 4, 5}, 4 ) ;" << '\n' << '\n' ;
    std::ptrdiff_t r11( 4 ) ;
    glas3::CP_tensor<glas3::dense_array<double>> t11( 0, s1, r11 ) ;
    glas3::randomize( t11, g ) ;

    std::cout << "tensor_train construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::tensor_train<glas3::dense_array<double>> t2(0, {5, 4, 3}, {4, 3}) ;" << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s2( {3, 4, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r2( {4, 3} ) ;
    glas3::tensor_train<glas3::dense_array<double>> t2(0, s2, r2) ;
    glas3::randomize( t2, g ) ;

    std::cout << "tensor_train construction (deep copy) from value, shape and rank" << '\n' ;
    std::cout << "glas3::tensor_train<glas3::dense_array<double>> t2(0, {5, 4, 3}, {2, 4}) ;" << '\n' ;
    glas3::dense_vector<std::ptrdiff_t> s22( {3, 4, 5} ) ;
    glas3::dense_vector<std::ptrdiff_t> r22( {2, 4} ) ;
    glas3::tensor_train<glas3::dense_array<double>> t22(0, s22, r22) ;
    glas3::randomize( t22, g ) ;

    std::cout << "entrywise_operation" << '\n' ;
    std::cout << "auto t111 = t1 + t11 + t1 ;" << '\n' ;
    auto t111 = t1 + t11 + t1 ;
    std::cout << "t111: " << t111 << '\n' ;
    glas3::randomize( t1, g ) ;
    glas3::randomize( t11, g ) ;
    auto t1_full = glas3::full( t1 ) ;
    auto t11_full = glas3::full( t11 ) ;
    auto t111_full = glas3::full( t111 ) ;
    auto t111_ref = t1_full + t11_full + t1_full ;
    if ( ! glas3::all( t111.shape() == t111_ref.shape() ) || glas3::norm_inf( t111_full - t111_ref ) > 1e-14 ) return 1 ;

    std::cout << "entrywise_operation" << '\n' ;
    std::cout << "auto t222 = t2 + t22 + t2 ;" << '\n' ;
    auto t222 = t2 + t22 + t2 ;
    std::cout << "t222: " << t222 << '\n' ;
    glas3::randomize( t2, g ) ;
    glas3::randomize( t22, g ) ;
    auto t2_full = glas3::full( t2 ) ;
    auto t22_full = glas3::full( t22 ) ;
    auto t222_full = glas3::full( t222 ) ;
    auto t222_ref = t2_full + t22_full + t2_full ;
    if ( ! glas3::all( t222.shape() == t222_ref.shape() ) || glas3::norm_inf( t222_full - t222_ref ) > 1e-14 ) return 1 ;

    std::cout << '\n' ;

	return 0;

}
