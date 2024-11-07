//  (C) Copyright Karl Meerbergen 2007. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#include "../../utils.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <random>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/timer.hpp>

typedef typename boost::uniform_real<double>                        distribution_type ;
typedef typename boost::lagged_fibonacci607                         engine_type ;
typedef boost::variate_generator< engine_type, distribution_type >  generator_type ;
std::random_device rd ;
generator_type g( engine_type( rd() ), distribution_type( -1, 1 ) ) ;

namespace ublas = boost::numeric::ublas ;

template <typename T>
struct do_test {
void operator() ( int N, std::ostream& f ) const {
  N = int( sqrt(double(N)) ) ;
  ublas::vector< T > v( N ), w( N ) ;
  ublas::matrix< T, ublas::column_major > a( N, N ) ;

  for ( int i = 0; i < N; ++i ) {
	  v(i) = g() ;
	  w(i) = g() ;
	  for ( int j = 0; j < N; ++j ) {
		  a(i, j) = g() ;
	  }
  }

  int const n_loop( (32/sizeof(T)) * (100000000 / (N*N) ) ) ;
  std::cout << "Doing " << n_loop << " loops of square matvec() for size " << v.size() << std::endl ;
  boost::timer timer ;
  for ( int i=0; i<n_loop; ++i ) {
    w.plus_assign( prod( trans(a), v ) ) ;
  }
  double time_elapsed = timer.elapsed() ;
  std::cout << "Sum : " << norm_2(w) << std::endl ;
  std::cout << "Elapsed time : " << time_elapsed << std::endl ;
  f << N << " " << ( time_elapsed * 1000000.0 ) / n_loop << std::endl ;
}
} ;
