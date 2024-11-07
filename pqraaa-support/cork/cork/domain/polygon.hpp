//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_polygon_hpp
#define cork_domain_polygon_hpp

#include <cork/domain/k_means.hpp>
#include <glas2/vector.hpp>
#include <cork/vector.hpp>
#include <cassert>
#include <limits>
#include <cmath>
#include <complex>
#include <cstring>
#include <type_traits>

namespace CORK { namespace domain {

  template <typename T>
  class polygon {
    public:
      typedef T                               value_type ;
      typedef decltype( std::abs(T()) )       real_type ;

      // dimension is 1 or 2
      typedef std::integral_constant<int,1>   dimension ;

    public:
      polygon( CORK::vector<value_type> const& points )
      : points_{points}
      , centre_( sum(points) / real_type(points_.size()) )
      {
        // Sort points following angles
      /*  glas2::vector<real_type> angles(points_.size()) ;
        for (int i=0; i<points_.size(); ++i) {
           angles(i) = std::acos( std::real( points_(i)-centre_ ) / std::abs( points(i)-centre_ ) ) ;
        }
        glas2::vector<int> order( points_.size() ) ; order = glas2::range(0,points_.size() ) ;
        glas2::sort( angles, order ) ;
        glas2::vector<value_type> points2( points_.size() ) ;
        points2 = points_(order) ;
        points_ = points2 ;*/
      }

    public:
     /**
      * Returns a glas vector containing +- n_points that discretized the domain.
      *
      * @param n_points   The amount of points approximately wanted.
      */
      auto discretize_coarse( int n_points ) const {
        return discretize( n_points ) ;
      } // discretize_coarse()


      auto /*const&*/ discretize( int n_points ) const {
        std::vector<value_type> points; points.reserve( n_points ) ;

        // Compute surface
        real_type surf = 0.0 ;
        for (int i=1; i<points_.size(); ++i) {
          surf += points_(i-1).real() * points_(i).imag() ;
          surf -= points_(i).real() * points_(i-1).imag() ;
        }
        surf += points_(points_.size()-1).real() * points_(0).imag() ;
        surf -= points_(0).real() * points_(points_.size()-1).imag() ;
        surf = 0.5 * std::abs( surf ) ;

        real_type cell_size = std::sqrt( surf / n_points ) ;
        assert(cell_size!=0.) ;

        // Discretize is small squares. Take the middle points
        real_type x_max = glas2::max( glas2::real(points_) ) ;
        real_type x_min = glas2::min( glas2::real(points_) ) ;
        real_type y_max = glas2::max( glas2::imag(points_) ) ;
        real_type y_min = glas2::min( glas2::imag(points_) ) ;
        for (real_type x=x_min; x<x_max; x+=cell_size) {
          for (real_type y=y_min; y<y_max; y+=cell_size) {
            value_type p( x, y ) ;
            //std::cout << p << " " << distance(p) << std::endl ;
            if (distance(p)==0.0) {
              points.push_back( p ) ;
            }
          }
        }
        glas2::vector<value_type> points_copy( points.size() ) ;
        std::copy( points.begin(), points.end(), points_copy.begin() ) ;
        return std::move(points_copy) ;
      } // discretize()



    public:
     /**
      * Returns the distance of a given point to the domain.
      *
      * @param point    A real or complex point of which the distance to the domain is calculated.
      */
      template <typename Point>
      auto distance( Point const& point ) const {
        // Check if point lies inside the polygon
        real_type distance = std::abs(point-points_(0)) ;
        int number = 0 ;
        for (int i=0; i<points_.size(); ++i) {
          value_type p1 ;
          value_type p2 ;
          if (i==0) {
            p1 = points_(points_.size()-1) ;
            p2 = points_(i) ;
          } else {
            p1 = points_(i) ;
            p2 = points_(i-1) ;
          }
          // Check if point lies between the two points.
          // If so, stop
          if (p1.imag()==p2.imag()) {
            if (p1.imag()==point.imag()) {
              if (point.real()<=p2.real() && point.real()>=p1.real()) {
                distance = 0.0 ;
                number = 1 ;
                break ;
              } else if (point.real()>=p2.real() && point.real()<=p1.real()) {
                distance = 0.0 ;
                number = 1 ;
                break ;
              }
            } else {
              if (point.real()<std::min(p1.real(),p2.real()) || point.real()>std::max(p1.real(),p2.real())) {
                real_type dist2 = std::min(std::abs(p1-point), std::abs(p2-point)) ;
                distance = std::min<real_type>(distance, dist2 ) ;
              } else {
                real_type dist2 = std::abs(p1.imag()-point.imag()) ;
                distance = std::min<real_type>(distance, dist2 ) ;
              }
            }
          } else {
            // Check intersection in horizontal positive direction
            real_type alpha = ( point.imag() - p2.imag() ) / (p1.imag()-p2.imag()) ;
            // use equality for 0, not for 1, otherwise intersection can be counted twice
            if (alpha>=0 && alpha<1.0) {
              real_type p = alpha*p1.real() + (1.-alpha)*p2.real() ;
              if (p==point.real()) {
                number = 1 ;
                break ;
              }
              if (p>point.real()) ++number ;
            }

            // Distance orthogonal to line segment.
            value_type q1( -p1.imag(), p1.real() ) ;
            value_type q2( -p2.imag(), p2.real() ) ;
            value_type v_alpha = p1-p2 ;
            value_type rhs = point - p2 ;
            real_type det = -std::norm(v_alpha) ;
            /*real_type*/ alpha = (-rhs.real()*v_alpha.real() - rhs.imag()*v_alpha.imag()) / det ;

            if (alpha<0.0) {
              real_type dist2 = std::abs( point - p2 ) ;
              distance = std::min( distance, dist2 ) ;
            } else if (alpha>1.0) {
              real_type dist2 = std::abs( point - p1 ) ;
              distance = std::min( distance, dist2 ) ;
            } else {
              real_type dist2 = std::abs( point - p1*alpha - p2*(1.-alpha) ) ;
              distance = std::min( distance, dist2 ) ;
            }
          }
        }
        if (number%2==1) distance = 0.0 ;
        return distance ;
      } // distance()

    private:
      // variables defining the domain
      CORK::vector<value_type> points_ ;
      value_type               centre_ ;
  } ;

} } // namespace CORK::domain

#endif
