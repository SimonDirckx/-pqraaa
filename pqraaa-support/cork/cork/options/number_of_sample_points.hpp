//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_NUMBER_OF_SAMPLE_POINTS_HPP
#define CORK_OPTIONS_NUMBER_OF_SAMPLE_POINTS_HPP

#include<limits>
#include<tuple>

namespace CORK { namespace options {

  class number_of_sample_points {
    public:
      inline number_of_sample_points()
      : value_(100000)
      {}

      template <typename ...Ts>
      number_of_sample_points( std::tuple<Ts...> const& options)
      : value_(100000)
      {}

      inline number_of_sample_points( int v )
      : value_(v)
      {}

      inline int value() const { return value_ ; }
      inline int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
