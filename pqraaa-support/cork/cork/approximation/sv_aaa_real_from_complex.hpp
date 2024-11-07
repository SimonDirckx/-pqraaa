//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_real_hpp
#define cork_approximation_sv_aaa_real_hpp

#include <cork/approximation/sv_aaa.hpp>
#include <cmath>

namespace CORK { namespace approximation {

  template <typename T, typename FunctionSequence, typename Domain>
  decltype(auto) SV_AAA_real( FunctionSequence const& fun_sequence, Domain const& domain ) {
    return SV_AAA_real( fun_sequence, domain, aaa_options< decltype(std::abs(T())) >() ) ;
  }

  template <typename T, typename FunctionSequence, typename Domain, typename Options>
  decltype(auto) SV_AAA_real( FunctionSequence const& fun_sequence, Domain const& domain, Options const& options ) {
    auto repr = SV_AAA<T>( fun_sequence, domain, options ) ;
    std::cout << repr.nodes() << "\n" ;

    // Duplicate all points by adding complex conjugate points
    // Count real points
    int n_real = 0 ;
    for (int i=0; i<repr.n(); ++i) {
      if (repr.nodes()(i).imag()==0.0) ++n_real ;
    }

    SV_AAA_triple< T, decltype(std::abs(T())) > repr_real( repr.coefficients().num_columns(), 2*repr.n()-n_real ) ;

    // First add the real points:
    int ii = 0 ;
    for (int i=0; i<repr.n(); ++i) {
      if (repr.nodes()(i).imag()==0.0) {
        repr_real.nodes_(ii) = repr.nodes()(i) ;
        repr_real.weights_(ii) = repr.weights()(i) ;
        assert( repr_real.weights_(ii).imag()==0.0 ) ;
        repr_real.coefficients_(ii, glas2::all()) = glas2::real( repr.coefficients()(i, glas2::all()) ) ;
        assert( norm_2( imag(repr.coefficients()(ii, glas2::all())) )==0.0 ) ;
        ++ii ;
      }
    }

    // Add the complex points:
    for (int i=0; i<repr.n(); ++i) {
      if (repr.nodes()(i).imag()!=0.0) {
        repr_real.nodes_(ii) = repr.nodes()(i) ;
        repr_real.weights_(ii) = repr.weights()(i) ;
        repr_real.coefficients_(ii, glas2::all()) = glas2::real( repr.coefficients()(i, glas2::all()) ) ;
        ++ii ;
        repr_real.nodes_(ii) = conj(repr.nodes()(i)) ;
        repr_real.weights_(ii) = conj(repr.weights()(i)) ;
        repr_real.coefficients_(ii, glas2::all()) = glas2::imag( repr.coefficients()(i, glas2::all()) ) ;
        ++ii ;
      }
    }

    return repr_real ;
  } // SV_AAA_real


} } // namespace CORK::approximation

#endif
