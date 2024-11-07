//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_options_initial_vector_hpp
#define cork_options_initial_vector_hpp

#include <cork/matrix.hpp>
#include <glas2/matrix.hpp>
#include <cork/krylov/user_defined_initial_vector.hpp>
#include <functional>
#include <iostream>

namespace CORK { namespace options {

  template <typename T>
  class initial_vector {
    public:
      typedef   T                 value_type;

      typedef   std::function<void( CORK::matrix<value_type> ) >    function_type;

    public:
      static initial_vector random() {
        return initial_vector( false, value_type(), [](CORK::matrix<value_type> v) {glas2::randomize(v);} ) ;
      }

      static initial_vector structured( value_type const& lambda, std::function<void( CORK::vector<value_type> ) > const& function ) {
        return initial_vector( true, lambda, [=]( CORK::matrix<value_type> v) {function(v(glas2::all(),0)) ;} ) ;
      }

      static initial_vector structured_random( value_type const& lambda ) {
        return initial_vector( true, lambda, [](CORK::matrix<value_type> v) {glas2::randomize(v(glas2::all(),0));} ) ;
      }

      static initial_vector user_defined( function_type const& function ) {
        return initial_vector( false, value_type(), function ) ;
      }

    public:
      initial_vector()
      : value_(false,value_type(),[](CORK::matrix<value_type> v) {glas2::randomize(v);})
      {}

      template <typename ...Ts>
      initial_vector( std::tuple<Ts...> const& options)
      : value_(false,value_type(), []( CORK::matrix<value_type> v ) {glas2::randomize(v);})
      {}

      initial_vector(bool structured, value_type const& lambda, function_type const& init_function) 
      : value_(structured, lambda, init_function)
      {}

    public:
      krylov::user_defined_initial_vector<value_type> const& value() const { return value_ ; }

  //  private:
 //     function_type default_function;

    private:
      krylov::user_defined_initial_vector<value_type> value_;
  };

} } // CORK::options

#endif
