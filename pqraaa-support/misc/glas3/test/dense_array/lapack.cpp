#include <glas3/bindings/array/dense_array/adaptor_dense_vector.hpp>
#include <glas3/bindings/array/dense_array/adaptor_dense_matrix.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>

#include <glas3/array/dense_array/algorithm/norm.hpp>
#include <glas3/array/dense_array/algorithm/ops.hpp>
#include <glas3/array/dense_array/algorithm/ttt.hpp>

#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>

int main() {
  using namespace boost::numeric::bindings ;

  //int n = 3 ;
  glas3::dense_vector<double> x({1, 2, 3}) ;
  glas3::dense_vector<int> ind({7, 8, 9}) ;
  glas3::dense_matrix<double> m({1, 2, 3, 0, 5, 9, -1, 13, 9}, {3, 3}) ;

  glas3::dense_vector<double> b = glas3::ttt(m, x, {1}, {0}) ;

  lapack::gesv( m, ind, b ) ;
  if ( glas3::norm_1( x - b ) > 1.e-10 * glas3::norm_1( x ) ) return 1 ;

//  m = 3.0 * glas2::eye( n, n ) ;
//  auto m_u = upper(m) ;
//  lapack::syev('N',m_u,x);

//  if (norm_2(x-glas2::constant(n,3.0))>1.e-10*norm_2(x)) return 2 ;

  return 0 ;
}
