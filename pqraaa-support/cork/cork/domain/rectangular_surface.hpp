//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_rectangular_surface_hpp
#define cork_domain_rectangular_surface_hpp

#include <glas2/vector.hpp>
#include <type_traits>
#include <cassert>
#include <complex>
#include <cmath>

namespace CORK { namespace domain {

  // Rectangle determined by 2 opposite corners
  template <typename T>
  class rectangular_surface {
    public:
      typedef T                             value_type ;
      typedef typename T::value_type        real_type ;
      typedef std::integral_constant<int,2> dimension ;

    public:
      rectangular_surface( T const& corner1, T const& corner4 )
      : corner_1_( corner1 )
      , corner_4_( corner4 )
      {
        assert( corner_1_.real()<corner_4_.real() ) ;
        assert( corner_1_.imag()<corner_4_.imag() ) ;
      }

    public:
      template <typename Points>
      Points discretize( Points points ) const {
        real_type ratio = std::sqrt(( corner_4_.real() - corner_1_.real() ) / ( corner_4_.imag() - corner_1_.imag() ) ) ;
        int n_imag = std::round( std::sqrt(points.size()) / ratio ) ;
        int n_real = (points.size()+n_imag-1) / n_imag ;
        real_type step_r = ( corner_4_.real() - corner_1_.real() ) ;
        if (n_real>1)
          step_r /= (n_real-1) ;
        real_type step_i = ( corner_4_.imag() - corner_1_.imag() ) ;
        if (n_imag>1)
          step_i /= (n_imag-1) ;
        int i_point = 0 ;
        value_type point = corner_1_ ;
        for (int i=0; i<n_imag; ++i) {
          for (int j=0; j<n_real; ++j, ++i_point) {
            if (i_point<points.size())
              points(i_point) = point ;
            point += step_r ;
          }
          point = value_type( corner_1_.real(), point.imag() + step_i ) ;
        }
        assert( i_point>=points.size() ) ;

        return points ;
      } // discretize()

      template <typename Points>
      void discretize( Points& points, real_type distance ) const {
        assert( distance>0. ) ;
        int n_points_r = std::floor( ( corner_4_.real() - corner_1_.real() ) / distance ) ;
        int n_points_i = std::floor( ( corner_4_.imag() - corner_1_.imag() ) / distance ) ;
        points.resize( n_points_r, n_points_i ) ;
        discretize( points ) ;
      } // discretize()

    public:
      template <typename Points>
      void discretize_coarse( Points points ) const {
        discretize( points ) ;
      }

      auto discretize( int n_points ) const {
        glas2::vector<value_type> points(n_points) ;
        discretize( points ) ;
        return std::move(points) ;
      }

    public:
      real_type distance( value_type const& point ) const {
        real_type hor_dist = std::max<real_type>( std::max(corner_1_.real() - point.real(), point.real()-corner_4_.real()), 0.0 ) ;
        real_type vert_dist = std::max<real_type>( std::max(corner_1_.imag() - point.imag(), point.imag()-corner_4_.imag()), 0.0 ) ;
        return hor_dist + vert_dist ;
      } // distance()

    private:
      value_type corner_1_ ;
      value_type corner_4_ ;
  } ;

} } // namespace CORK::domain

#endif
