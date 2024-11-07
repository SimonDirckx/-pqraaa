//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_half_disk_hpp
#define cork_domain_half_disk_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <cmath>
#include <complex>
#include <type_traits>

namespace CORK { namespace domain {

  // Rectangle determined by 2 opposite corners
  template <typename T>
  class half_disk {
    public:
      typedef T                             value_type ;
      typedef typename T::value_type        real_type ;
      typedef std::integral_constant<int,2> dimension ;

    public:
      // If right_end_point lies on the left, then the disk is oriented towards below.
      // i.e., the right_end_point is the first point, when walking around the disk counter clockwise, starting with angle zero.
      half_disk( value_type const& center, value_type const& right_end_point )
      : right_end_point_( right_end_point )
      , center_( center )
      , radius_( std::abs(right_end_point_-center_) )
      , axis_dir_( (right_end_point_- center_)/std::abs(right_end_point_- center_) )
      {
        assert( right_end_point_!=center_ ) ;
      } 

/*    public:
      bool contains( value_type const& point ) const {
        value_type dir_p = (point-center_)/std::abs(point-center_) ;
        return ( std::abs(point-center_) <=radius_ )
          &&   ( std::imag( dir_p / axis_dir_ ) > 0. )
          ;
      } // contains()

      template <typename Real>
      bool contains( value_type const& point, Real const& tolerance ) const {
        value_type dir_p = (point-center_)/std::abs(point-center_) ;
        return ( std::abs(point-center_) <=radius_+tolerance*radius_ )
          &&   ( std::imag( dir_p / axis_dir_ ) > -tolerance )
          ;
      } // contains()

      real_type surface() const {
        return std::atan(1.0)*2.0 * radius_*radius_ ;
      } // surface()
*/
    private:
      template <typename Points>
      void discretize( Points points ) const {
        assert(points.size()>1) ;
        real_type const pi = 4.*std::atan(1.0) ;
        real_type surf = pi * radius_ * radius_ ;
        real_type element = surf / points.size() ;
        real_type radial_step = radius_ * std::sqrt( element/surf ) ;

        int n_curves = std::floor(radius_ / radial_step) ;
        // total length of all curves:
        real_type total_length = pi * radial_step*(n_curves+1)*n_curves*0.5 ;

        // Actual step on an arc
        real_type step = total_length / points.size() ;

        // Number of points on longest half circle:
        int n_point_max = std::floor( pi * radius_ / step ) ;

        int i_point = 0 ;
        points(i_point) = center_ ;
        ++i_point ;
        real_type radius = radial_step ;
        for (int i=0; i<n_curves; ++i, radius+=radial_step) {
          // Precise division over arc:
          int m_points = std::round( n_point_max / radius_ * radius ) ;
          real_type angle_step = pi / (m_points-1) ;
          assert( angle_step>0. ) ;
          if (i==n_curves-1) angle_step = pi / (points.size()-i_point-1) ;
          for (real_type angle=0*0.5*angle_step; angle<pi+0.5*angle_step; angle+=angle_step ) {
            if (i_point<points.size())
              points( i_point ) = center_ + axis_dir_ * radius * value_type( std::cos(angle), std::sin(angle) ) ;
            ++i_point ;
          }
        }
      } // discretize()

    public:
      glas2::vector<value_type> discretize( int n_points ) const {
        glas2::vector<value_type> points( n_points ) ;
        discretize( points ) ;
        return std::move(points) ;
      } // discretize()

      template <typename Points>
      void discretize_coarse( Points points ) const {
        discretize( points ) ;
      } // discretize()


      real_type distance( value_type const& point ) const {
        return 0. ;
      }

    private:
      value_type right_end_point_ ;
      value_type center_ ;
      real_type  radius_ ;
      value_type axis_dir_ ;
  } ;

} } // namespace CORK::domain

#endif
