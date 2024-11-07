//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)
#ifndef cork_options_stop_criterion_eigen_vec_needed_hpp
#define cork_options_stop_criterion_eigen_vec_needed_hpp

#include <cork/options/stop_criterion.hpp>
#include <tuple>

// should be obsolete
namespace CORK { namespace options {

  template <typename T>
  class stop_criterion_eigen_vec_needed {

    public:
      stop_criterion_eigen_vec_needed()
      : value_(false)
      {}

      template <typename ...Ts>
      stop_criterion_eigen_vec_needed( std::tuple<Ts...> const& options )
      : value_(options::value_of<options::stop_criterion<T>>(options).eig_vec_needed())
      {}

      stop_criterion_eigen_vec_needed( bool value )
      :value_(value)
      {}

    public:
      bool  value() const { return value_ ; }
      bool& value()       { return value_ ; }

    private:
      bool      value_;

  } ; // class stop_criterion_eigen_vec_needed

} } // CORK::options

#endif
