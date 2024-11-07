//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_SAMPLE_FUNCTIONS_HPP
#define CORK_OPTIONS_SAMPLE_FUNCTIONS_HPP

#include<tuple>
#include<cork/basis/set_of_functions.hpp>
#include<cork/vector.hpp>
#include<functional>

namespace CORK { namespace options {

  template <typename T>
  class sample_functions {
    public:
      typedef std::function<void(T,CORK::vector<T>)> function_type ;

    public:
      sample_functions()
      : value_( []( T, CORK::vector<T> v ) { assert(false) ; }, 0 )
      {}

      template <typename ...Ts>
      sample_functions( std::tuple<Ts...> const& options)
      : value_( []( T, CORK::vector<T> v ) { assert(false) ; }, 0 )
      {}

      sample_functions( function_type const& f, int n_functions )
      : value_( f, n_functions )
      {}

      auto value() const { return value_ ; }

    private:
      basis::set_of_functions< T, function_type > value_ ;
  } ;

} } // CORK::options

#endif
