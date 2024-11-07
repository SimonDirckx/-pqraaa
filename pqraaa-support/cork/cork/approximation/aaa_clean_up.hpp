//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_clean_up_hpp
#define cork_approximation_sv_aaa_clean_up_hpp

#include <cork/approximation/aaa_options.hpp>
#include <cork/approximation/sv_aaa_approximation.hpp>
#include <cork/approximation/sv_aaa_2_partial_fractions.hpp>
#include <glas2/bindings/vector.hpp>

namespace CORK { namespace approximation {

  template <typename AAAApproximation>
  decltype(auto) sv_aaa_clean_up( AAAApproximation const& aaa, aaa_options const& options ) {
    // Determine poles
    auto pf = sv_aaa_2_partial_fractions( aaa ) ;

    std::vector< value_type > supports( aaa.n() ) ;
    std::vector< real_type > supports_transformed( aaa.n() ) ;
    glas2::bind_vector(supports) = aaa.nodes() ;

    // Find small residues
    for (int i=0; i<pf.n(); ++i) {
      if (norm_inf(pf.coefficients()(i, glas2::all()))<options.tolerance) {
        // remove nearest support point
        glas2::bind_vector(supports_transformed) = glas2::abs( glas2::bind_vector(supports) - pf.nodes()(i) ) ;
        auto it = std::min_element( supports_transformed.begin(), supports_transformedupport.end() ) ;
        supports.remove( supports.begin()+(it-supports_transformed.begin()) ) ;
        supports_transformed.remove( it ) ;
      }
    }

    glas2::shared_vector< typename AAAApproximation::value_type > points( supports.size() ) ;
    points = glas2::bind_vector(supports) ;
  } // sv_aaa_clean_up()

} } // namespace CORK::approximation

#endif
