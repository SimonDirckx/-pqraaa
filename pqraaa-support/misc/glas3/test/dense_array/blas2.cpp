#include <glas3/bindings/array/dense_array/adaptor_dense_vector.hpp>
#include <glas3/bindings/array/dense_array/adaptor_dense_matrix.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>

#include <glas3/array/dense_array/algorithm/norm.hpp>
#include <glas3/array/dense_array/algorithm/ops.hpp>
#include <glas3/array/dense_array/algorithm/ttt.hpp>

#include <boost/numeric/bindings/blas/level2.hpp>

int main() {

  using namespace boost::numeric::bindings ;

  glas3::dense_vector<double> v({1, 2, 3}), w({4, 5}) ;
  glas3::dense_vector<double> x = w ;
  glas3::dense_matrix<double> m({1, 2, 3, 4, 5, 6}, {2, 3}) ;

  blas::gemv( 2.0, m, v, 1.0, w ) ;
  if ( glas3::norm_1( w - x - glas3::dense_vector<double>(2.0, 2) * glas3::ttt(m, v, {1}, {0}) ) > 1.e-12 ) return 1 ;

  return 0 ;

}
