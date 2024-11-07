//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_k_means_hpp
#define cork_domain_k_means_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <limits>
#include <cmath>
#include <cstring>
#include <type_traits>

namespace CORK { namespace domain {

  template <typename Points, typename Clusters>
  auto k_means( Points const& points, Clusters clusters, int max_it ) {
    typedef typename Points::value_type        value_type ;
    typedef decltype( std::abs(value_type()) ) real_type ;

    assert( clusters.size()<=points.size() ) ;
    glas2::vector< int > cluster_for_point( points.size() ) ;
    glas2::vector< value_type > old_clusters( clusters.size() ) ;
    glas2::vector< int > cluster_sizes( clusters.size() ) ;
    fill( cluster_for_point, -1 ) ;

    for (int i=0; i<clusters.size(); ++i) {
      glas2::range r_i( points.size()/clusters.size()*i, std::min(points.size()/clusters.size()*(i+1),points.size()) ) ;
      old_clusters(i) = glas2::sum( points( r_i ) ) / real_type(r_i.size()) ;
      old_clusters(i) = points( r_i(0) ) ;
    }

    bool done = false ;
    for (int i=0; i<max_it && !done; ++i) {
      fill(clusters, 0.0) ;
      fill(cluster_sizes, 0) ;
      done = true ;
      for (int j=0; j<points.size(); ++j) {
        int cluster = glas2::max_ind( - glas2::abs_squared( points(j) - old_clusters ) ) ;
        if (cluster!=cluster_for_point(j)) done = false ;
        cluster_for_point(j) = cluster ;
        clusters( cluster ) += points(j) ;
        cluster_sizes( cluster ) ++ ;
      }

      clusters /= cluster_sizes ;
      old_clusters = clusters ;
    } // for(i)

    // Take points as cluster points.
    for (int i=0; i<clusters.size(); ++i) {
      clusters(i) = points( glas2::max_ind( - glas2::abs_squared(clusters(i)-points) ) ) ; 
    }
  } // k_means()

} } // namespace CORK::domain

#endif
