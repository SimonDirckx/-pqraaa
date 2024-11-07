//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_point_cloud_hpp
#define cork_domain_point_cloud_hpp

#include <glas2/vector.hpp>
#include <cork/containers/vector.hpp>
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

    public:
      point_cloud( CORK::vector<value_type> points )
        : points_{points}
      {
        glas2::sort_comp(points_, realcomp); 
      }

    public:
      static bool realcomp(value_type a, value_type b) {
        if (std::real(a) == std::real(b)) return std::imag(a) < std::imag(b);
        return std::real(a) < std::real(b);
      } 

    public:
     /**
      * Returns a glas vector containing +- n_points that discretized the domain.
      *
      * @param n_points   The amount of points approximately wanted.
      */
      template <typename U>
      typename std::enable_if< std::is_integral<U>::value, glas2::shared_vector<value_type> >::type
      discretize( U n_points ) const {
        if (n_points <= points_.size() ) { 
          glas2::shared_vector<value_type> points(n_points);
          int step = (points_.size()-1)/(n_points-1);
          int startval = (points_.size()-(n_points-1)*step-1)/2;
          for (int i=startval, j=0; i<points_.size(); i+=step, j++) {
            points(j) = points_(i);
          }
          return points;
        } else {
          return points_;
        }
      } // discretize()
    
     /**
      * Fills a glas densevector with discretized points and returns the amount of discretized points.
      *
      * The amount of points discretized can differ from the size of the vector, but should be maximum
      * this amount.
      *
      * @param points   A glas2 dense vector which will be filled with discretized points.
      */
      template <typename Points>
      typename std::enable_if< glas2::is< glas2::DenseVector, typename Points::size_type>::value, Points>::type 
      discretize( Points points ) const {
        for (typename Points::size_type i=0; i<points.size(); ++i) {
          points(i) = points_(i);
        } 
        return points.size() ;
      } // discretize()

    public:
     /**
      * Returns the distance of a given point to the domain.
      *
      * @param point    A real or complex point of which the distance to the domain is calculated.
      */
      template <typename Point>
      auto distance( Point const& point ) const {
        real_type min_dist=std::numeric_limits<value_type>::max();
        for (int i=0; i<points_.size(); i++) {
          if (std::abs(points_(i)-point) < min_dist) {
            min_dist = std::abs(points_(i)-point);
          }
        }
        return min_dist;
      } // distance()

    private:
      // variables defining the domain
      const glas2::shared_vector<value_type> points_;

  } ;

} } // namespace CORK::domain

#endif
