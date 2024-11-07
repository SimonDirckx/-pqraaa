//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_ellipse_hpp
#define cork_domain_ellipse_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <cmath>
#include <complex>

namespace CORK { namespace domain {

  // Ellipse, defined by two focal points and the length of thelength of the semi major axis (the radius of the axis going through the two focal points)
  template <typename T>
  class ellipse {
    public:
      typedef T                              value_type ;
      typedef decltype(std::abs(value_type())) real_type ;
      typedef std::integral_constant<int,1> dimension ;

    public:
      ellipse( value_type const& focus1, value_type const& focus2, real_type const& semi_major_axis )
      : semi_major_axis_( semi_major_axis )
      , focus_1_( focus1 )
      , focus_2_( focus2 )
      {
        semi_minor_axis_ = std::pow(semi_major_axis_,2) - 0.25*std::norm(focus_2_-focus_1_) ;
        assert( semi_minor_axis_ > 0.0 ) ;
        semi_minor_axis_ = std::sqrt( semi_minor_axis_ ) ;
      } 

    public:
      bool contains( value_type const& point ) const {
        return std::abs(point-focus_1_) + std::abs(point-focus_2_)<=semi_major_axis_*2.0 ;
      } // contains()

      template <typename Real>
      bool contains( value_type const& point, Real const& tolerance ) const {
        return std::abs(point-focus_1_) + std::abs(point-focus_2_)<=(semi_major_axis_*2.0+tolerance) ;
      } // contains()

    private:
      /*template <typename Points>
      void discretize_deterministic( Points points ) const {
        assert(points.size()>1) ;
        real_type const pi = 4.*std::atan(1.0) ;
        real_type surf = pi * semi_major_axis_ * semi_minor_axis_ ;
        real_type element = surf / points.size() ;
        real_type radial_step = semi_major_axis_ * std::sqrt( element/surf ) ;
        real_type small_radial_step = semi_minor_axis_ * std::sqrt( element/surf ) ;

        int n_curves = std::floor(semi_major_axis_ / radial_step) ;
        // total length of all curves:
        real_type h = std::pow( (semi_major_axis_-semi_minor_axis_)/(semi_minor_axis_+semi_major_axis_), 2 ) ;
        real_type total_length = pi * (semi_minor_axis_+semi_major_axis_) * (n_curves+1) * (1.0 + 3*h/(10.+std::sqrt(4.-3.*h))) ;

        // Actual step on an arc
        real_type step = total_length / points.size() ;

        // Number of points on longest ellipse:
        int n_point_max = 2.*total_length / step / (n_curves+1) ;

        value_type center = 0.5*(focus_1_+focus_2_) ;
        value_type direction = focus_2_-focus_1_ ;
        if (direction==0.0)
          direction = 1.0 ;
        else
          direction /= std::abs( direction ) ;

        int i_point = 0 ;
        points(i_point) = focus_1_ ;
        ++i_point ;
        points(i_point) = focus_2_ ;
        ++i_point ;
        real_type radius = radial_step ;
        real_type small_radius = small_radial_step ;
        for (int i=0; i<n_curves; ++i, radius+=radial_step, small_radius+=small_radial_step) {
          // Precise division over arc:
          int m_points = std::round( n_point_max / semi_major_axis_ * radius ) ;
          real_type angle_step = 2*pi / (m_points-1) ;
          assert( angle_step>0. ) ;
          if (i==n_curves-1) angle_step = 2*pi / (points.size()-i_point-1) ;
          for (real_type angle=0.5*angle_step; angle<2*pi+0.5*angle_step; angle+=angle_step ) {
            if (i_point<points.size())
              points( i_point ) = center + direction * value_type( radius*std::cos(angle), small_radius*std::sin(angle) ) ;
            ++i_point ;
          }
        }
      } // discretize_deterministic()*/

    public:
      auto discretize_coarse( int n_points ) const {
        return discretize( n_points ) ;
      }

      auto discretize( int n_points ) const {
        glas2::vector<value_type> points( n_points ) ;
        discretize_impl( points ) ;
        return std::move(points) ;
      } // discretize()

    private:
      template <typename Points>
      void discretize_impl( Points points ) const {
        real_type const pi = 4.*std::atan(1.0) ;
        real_type surf = pi * semi_major_axis_ * semi_minor_axis_ ;
        real_type element = surf / points.size() ;
        real_type step = std::sqrt( element ) ;

        value_type center = 0.5*(focus_1_+focus_2_) ;
        value_type direction = focus_2_-focus_1_ ;
        if (direction==0.0) direction=1.0;
        direction /= std::abs(direction) ;
        value_type orto_direction = value_type( direction.imag(), direction.real() ) ;

        value_type point ;

        int i_pt = 0 ;
        for (real_type vert=0.5*step; vert<semi_minor_axis_; vert+=step) {
          //for ( real_type hor=semi_major_axis_*std::cos(std::asin(vert/semi_minor_axis_)); hor>=0.0; hor-=step) {
          real_type max_hor = semi_major_axis_*std::sqrt(semi_minor_axis_*semi_minor_axis_-vert*vert)/semi_minor_axis_;
          for ( real_type hor=0.0; hor<=max_hor; hor+=step) {
            point = center + direction*hor + orto_direction*vert ;
            if (i_pt<points.size()) points(i_pt) = point ; ++i_pt ;
            point = center - direction*hor + orto_direction*vert ;
            if (i_pt<points.size()) points(i_pt) = point ; ++i_pt ;
            point = center - direction*hor - orto_direction*vert ;
            if (i_pt<points.size()) points(i_pt) = point ; ++i_pt ;
            point = center + direction*hor - orto_direction*vert ;
            if (i_pt<points.size()) points(i_pt) = point ; ++i_pt ;
          }
        }
        assert( i_pt>=points.size() ) ;
      } // discretize()

    public:
      real_type distance( value_type const& point ) const {
        value_type center = 0.5*(focus_1_+focus_2_) ;
        value_type direction = focus_2_-focus_1_ ;
        if (direction==0.0) direction=1.0;
        direction /= std::abs(direction) ;
        value_type rel_point = point - center ;
        real_type cos_d = std::real(direction) ;
        real_type sin_d = std::imag(direction) ;
        real_type abs_p = std::abs(rel_point) ;
        real_type cos_p = std::real(rel_point) / abs_p ;
        real_type sin_p = std::imag(rel_point) / abs_p ;
        real_type cos_t = cos_p*cos_d + sin_p*sin_d ; cos_t = cos_t * cos_t ;
        real_type sin_t = 1. - cos_t ;
        real_type rel_point_x = abs_p * abs_p * cos_t ;
        real_type rel_point_y = abs_p * abs_p * sin_t ;
        return std::max( 0.0, rel_point_x / (semi_major_axis_*semi_major_axis_) + rel_point_y /(semi_minor_axis_*semi_minor_axis_) -1.0 ) ;
      }

    private:
      real_type  semi_major_axis_ ;
      value_type focus_1_ ;
      value_type focus_2_ ;
      real_type  semi_minor_axis_ ;
  } ;

/*
  template <typename T>
  auto ellipse( T const& focus1, T const& focus2, decltype(std::abs(T())) const& semi_major_axis ) {
    return ellipse_type<T>( focus1, focus2, semi_major_axis ) ;
  }
*/
} } // namespace CORK::domain

#endif
