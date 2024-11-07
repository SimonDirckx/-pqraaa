#ifndef cork_utility_norm_est_hpp
#define cork_utility_norm_est_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <complex>
#include <limits>
#include <utility>

namespace CORK {

  //
  // Compute an estimate of the two norm by matrix vector products using the function
  //   F( x, y )
  //
  // Example:
  //   auto nrm = norm_est<T>( n, [](auto const& x, auto y ) { ... } ;
  //
  template <typename T, typename F>
  decltype (auto) norm_est( int n, F const& f ) {
    glas2::shared_vector< T > v( n ) ;
    glas2::shared_vector< T > w( n ) ;
    randomize( v ) ; v /= norm_2(v) ;
    f( v, w ) ;
    auto mu = norm_2( w ) * sqrt( decltype(std::abs(T()))(n) ) ;
    return std::make_pair( mu / 3.0, mu * 100.0 ) ;
  } // norm_est()

} // namespace CORK

#endif
