//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_rectangle_hpp
#define cork_domain_rectangle_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <complex>
#include <cmath>
#include <type_traits>

namespace CORK { namespace domain {

  // Rectangle determined by 2 opposite corners
  template <typename T>
  class rectangle {
    public:
      typedef T                             value_type ;
      typedef typename T::value_type        real_type ;
      typedef std::integral_constant<int,1> dimension ;

    public:
      rectangle( T const& corner1, T const& corner4 )
      : corner_1_( corner1 )
      , corner_4_( corner4 )
      {
        assert( corner_1_.real()<corner_4_.real() ) ;
        assert( corner_1_.imag()<corner_4_.imag() ) ;
      }
/*
    private:
      template <typename Points>
      Points discretize( Points points ) const {
        real_type length = corner_4_.real()-corner_1_.real() + corner_4_.imag()-corner_1_.imag() ;

        // Number for first part:
        int first_halve = points.size() / 2 ;
        int first_quarter = round( first_halve * real(corner_4_-corner_1_) / length ) ;
        int i_point = 0 ;
        real_type step = real(corner_4_-corner_1_) / first_quarter ;
        value_type point = corner_1_ ;
        for ( ; i_point<first_quarter; ++i_point ) {
          point += step ;
          points(i_point) = point ;
        }
        point = value_type( corner_4_.real(), corner_1_.imag() ) ;
        step = imag(corner_4_-corner_1_) / (first_halve-i_point) ;
        for ( ; i_point<first_halve; ++i_point ) {
          point = value_type( point.real(), imag(point) + step ) ;
          points(i_point) = point ;
        }

        first_quarter += first_halve ;
        point = corner_4_ ;
        step = real(corner_4_-corner_1_) / (first_quarter-i_point) ;
        for ( ; i_point<first_quarter; ++i_point ) {
          point -= step ;
          points(i_point) = point ;
        }

        point = value_type( corner_1_.real(), corner_4_.imag() ) ;
        step = imag(corner_4_-corner_1_) / (points.size()-i_point) ;
        for ( ; i_point<points.size(); ++i_point ) {
          point = value_type( point.real(), point.imag()-step ) ;
          points(i_point) = point ;
        }
        assert(i_point==points.size()) ;

        return points ;
      } // discretize()
      */

    public:
      auto discretize_coarse( int n_points ) const {
        return discretize( n_points ) ;
      }

      auto discretize( int n_points ) const {
        glas2::vector<value_type> points( n_points ) ;
        discretize_impl( points ) ;
        return std::move(points) ;
      } // discretize()

    public:
      real_type distance( value_type const& point ) const {
        real_type hor_dist = std::max(corner_1_.real() - point.real(), point.real()-corner_4_.real()) ;
        real_type vert_dist = std::max(corner_1_.imag() - point.imag(), point.imag()-corner_4_.imag()) ;
        if (hor_dist<0.0 && vert_dist<0.0) return std::abs( std::max( hor_dist, vert_dist ) ) ;
        if (hor_dist<=0.0) return vert_dist ;
        if (vert_dist<=0.0) return hor_dist ;
        return hor_dist + vert_dist ;
      } // distance()

    private:
      value_type corner_1_ ;
      value_type corner_4_ ;
  } ;

} } // namespace CORK::domain

#endif
