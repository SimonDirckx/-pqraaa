//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../../utils.hpp"
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas ;

template <typename T>
struct do_test {
void operator() ( int N, std::ostream& f ) const {
  N = int( sqrt(double(N)) ) ;
  ublas::vector< T > v( N ), w( N ) ;
  ublas::matrix< T, ublas::column_major > a( N, N ) ;

  glas2::seed< decltype( std::abs(T()) ) >  seed ;
  glas2::randomize( glas2::bind_vector(v), seed ); glas2::randomize( glas2::bind_vector(w), seed ) ;
  for (int i=0; i<N; ++i) glas2::randomize( glas2::bind_matrix(a), seed );

  int const n_loop( (32/sizeof(T)) * (100000000 / (N*N) ) ) ;
  std::cout << "Doing " << n_loop << " loops of square matvec() for size " << v.size() << std::endl ;
  boost::timer::cpu_timer timer ;
  for ( int i=0; i<n_loop; ++i ) {
    w.plus_assign( prod( trans(a), v ) ) ;
  }
  double time_elapsed = timer.elapsed().wall ;
  std::cout << "Sum : " << norm_2(w) << std::endl ;
  std::cout << "Elapsed time : " << time_elapsed << std::endl ;
  f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
}
} ;
