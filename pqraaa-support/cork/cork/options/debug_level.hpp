//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_DEBUG_LEVEL_HPP
#define CORK_OPTIONS_DEBUG_LEVEL_HPP

#include<tuple>

namespace CORK { namespace options {

  class debug_level {
    public:
      debug_level()
      : value_(0)
      {}

      debug_level(int level)
      : value_(level)
      {}

      template <typename ...Ts>
      debug_level( std::tuple<Ts...> const& options)
      : value_(0)
      {}

      int value() const { return value_ ; }
      int& value() { return value_ ; }

    private:
      int value_ ;
  } ;

} } // CORK::options

#endif
