//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_line_segment_hpp
#define cork_domain_line_segment_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <limits>
#include <cmath>
#include <complex>
#include <algorithm>
#include <cstring>
#include <type_traits>

namespace CORK { namespace domain {

  template <typename T>
  class line_segment {
    public:
      typedef T                               value_type ;
      typedef decltype( std::abs(T()) )       real_type ;

    public:
      line_segment( T const& a, T const& b )
      : a_(a)
      , b_(b)
      {} // line_segment() constructor

    public:
     /**
      * Returns a glas vector containing +- n_points that discretized the domain.
      *
      * @param n_points   The amount of points approximately wanted.
      */
      auto discretize( int n_points ) const {
        glas2::vector<value_type> points(n_points);
        points = glas2::linspace( a_, b_, n_points );
        return std::move(points) ;
      } // discretize()

      auto discretize_coarse( int n_points ) const {
        return discretize( n_points ) ;
      } // discretize()

    public:
/*      //// NOT YET OBSOLETE ==>
      template <typename Point>
      bool contains( Point const& point ) const {
        return contains( point, std::numeric_limits< decltype(std::abs(value_type())) >::epsilon() ) ;
      } // contains()

      template <typename Real, typename Point>
      bool contains( Point const& point, Real const& tolerance ) const {
        auto z = (point-a_)/(b_-a_);
        if (std::real(z) >= 0 && std::real(z) <= 1) return std::abs(std::imag(z)*(b_-a_)) <= tolerance;
        return std::min(std::abs(point-a_), std::abs(point-b_)) <= tolerance; 
      } // contains()
      //// <== NOT YET OBSOLETE
      */

     /**
      * Returns the distance to the domain.
      *
      * @param point  The point of which the distance is computed.
      */
      template <typename Point>
      auto distance( Point const& point ) const {
        typedef std::common_type_t<Point, value_type> ctype ;
        ctype a = a_;
        ctype b = b_;
        ctype p = point;
        // linear transformation of the space which transforms the line [a_ b_]
        // to the line [0+0i, 1+0i] which makes it easier to chose the distance formula 
        ctype z = (p-a)/(b-a);
        // Real value smaller than 0 => distance to point a;
        if (std::real(z) <= 0) { return std::abs(p-a);}
        // Real value greater than 1 => distance to point b;
        if (std::real(z) >= 1) { return std::abs(p-b);}
        // else distance to line
        return std::abs(std::imag(z)*(b-a));
      }

      value_type const& a() const {  return a_ ; }
      value_type const& b() const {  return b_ ; }

    private:
      // variables defining the domain
      value_type a_;
      value_type b_;

  } ;

/*
  template <typename T>
  auto line_segment( T const& a, T const& b ) {
    return line_segment_type<T>( a, b ) ;
  }
*/
} } // namespace CORK::domain

#endif
