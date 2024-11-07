//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_point_cloud_hpp
#define cork_domain_point_cloud_hpp

#include <cork/domain/k_means.hpp>
#include <glas2/vector.hpp>
#include <cork/vector.hpp>
#include <cassert>
#include <limits>
#include <cmath>
#include <cstring>
#include <type_traits>

namespace CORK { namespace domain {

  template <typename T>
  class point_cloud {
    public:
      typedef T                               value_type ;
      typedef decltype( std::abs(T()) )       real_type ;

      // dimension is 1 or 2
      typedef std::integral_constant<int,1>   dimension ;

    public:
      point_cloud( CORK::vector<value_type> const& points )
        : points_{points}
      {
        glas2::sort_comp(points_, realcomp); 
      }

    public:
      static bool realcomp(value_type a, value_type b) {
        if (glas2::real(a) == glas2::real(b)) return glas2::imag(a) < glas2::imag(b);
        return glas2::real(a) < glas2::real(b);
      } 

    public:
     /**
      * Returns a glas vector containing +- n_points that discretized the domain.
      *
      * @param n_points   The amount of points approximately wanted.
      */
      auto discretize_coarse( int n_points ) const {
        assert( n_points<=points_.size() ) ;
        glas2::shared_vector< value_type > points( n_points ) ;
        k_means( points_, points, std::max<int>(10,points.size()/2) ) ;
        return points ;
      } // discretize()


      auto /*const&*/ discretize( int n_points ) const {
        return points_ ;
      } // discretize()



     /**
      * Fills a glas densevector with discretized points and returns the amount of discretized points.
      *
      * The amount of points discretized can differ from the size of the vector, but should be maximum
      * this amount.
      *
      * @param points   A glas2 dense vector which will be filled with discretized points.
      * /
      template <typename Points>
      typename std::enable_if< glas2::is< glas2::DenseVector, typename Points::size_type>::value, Points>::type 
      discretize( Points points ) const {
        for (typename Points::size_type i=0; i<points.size(); ++i) {
          points(i) = 0;
        } 
        return points.size() ;
      } // discretize()
      */

    public:
     /**
      * Returns the distance of a given point to the domain.
      *
      * @param point    A real or complex point of which the distance to the domain is calculated.
      */
      template <typename Point>
      auto distance( Point const& point ) const {
        return std::sqrt( - glas2::max( - glas2::abs_squared(point-points_) ) ) ;
      } // distance()

    private:
      // variables defining the domain
      const CORK::vector<value_type> points_;

  } ;
/*
  template <typename T>
  auto point_cloud( CORK::vector<T> const& points ) {
    return point_cloud_type<T>( points ) ;
  }*/

} } // namespace CORK::domain

#endif
