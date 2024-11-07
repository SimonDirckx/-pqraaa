//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_TEST_APPROXIMATION_ERROR_HPP
#define CORK_OPTIONS_TEST_APPROXIMATION_ERROR_HPP

#include<tuple>

namespace CORK { namespace options {

  class test_approximation_error {
    public:
      test_approximation_error()
      : value_( false )
      {}

      template <typename ...Ts>
      test_approximation_error( std::tuple<Ts...> const& options)
      : value_( false )
      {}

      test_approximation_error( bool v )
      : value_(v)
      {}

      bool value() const { return value_ ; }
      bool& value() { return value_ ; }

    private:
      bool value_ ;
  } ;

} } // CORK::options

#endif
