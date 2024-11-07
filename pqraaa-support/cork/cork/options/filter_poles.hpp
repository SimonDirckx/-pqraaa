//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_FILTER_POLES_HPP
#define CORK_OPTIONS_FILTER_POLES_HPP

#include<tuple>
#include<cork/vector.hpp>
#include<functional>

namespace CORK { namespace options {

  template <typename T>
  class filter_poles {
    public:
      typedef std::function<int(CORK::vector<T>)> function_type ;

    public:
      filter_poles()
      : value_( []( CORK::vector<T> v ) { assert( false ) ; return v.size() ; } )
      {}

      template <typename ...Ts>
      filter_poles( std::tuple<Ts...> const& options)
      : value_( []( CORK::vector<T> v ) {assert( false ) ;  return v.size() ; } )
      {}

      filter_poles( function_type const& f )
      : value_( f )
      {}

      function_type value() const { return value_ ; }
      function_type& value() { return value_ ; }

    private:
      function_type value_ ;
  } ;

} } // CORK::options

#endif
