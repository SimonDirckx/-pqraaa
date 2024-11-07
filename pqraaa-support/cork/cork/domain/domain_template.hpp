//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_domain_domain_template_hpp
#define cork_domain_domain_template_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <limits>
#include <cmath>
#include <cstring>
#include <type_traits>

namespace CORK { namespace domain {

  template <typename T>
  class domain_template {
    public:
      typedef T                               value_type ;
      typedef decltype( std::abs(T()) )       real_type ;

      // dimension is 1 or 2
      typedef std::integral_constant<int,1>   dimension ;

    public:
      domain_template(  )
      {}

    public:
     /**
      * Returns a glas vector containing +- n_points that discretized the domain.
      *
      * @param n_points   The amount of points approximately wanted.
      */
      template<typename U, typename std::enable_if<std::is_integral<U>::value>::type* = nullptr>
      glas2::vector<value_type> discretize( U n_points ) const {
        // Can do something smarter here where more or less than n_points are used.
        glas2::vector<value_type> points( n_points ) ;
        discretize( points ) ;
        return std::move(points) ;
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
        for (::size_type i=0; i<points.size(); ++i) {
          points(i) = 0;
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
        return 0;
      } // distance()

    private:
      // variables defining the domain

  } ;

} } // namespace CORK::domain

#endif
