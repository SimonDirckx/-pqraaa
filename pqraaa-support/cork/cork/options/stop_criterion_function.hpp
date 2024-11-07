//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)
#ifndef CORK_OPTIONS_STOP_CRITERION_FUNCTION_HPP
#define CORK_OPTIONS_STOP_CRITERION_FUNCTION_HPP

#include <cork/containers/vector.hpp>
#include <functional>
#include <limits>
#include <tuple>
#include <complex>

namespace CORK { namespace options {

  template <typename T>
  struct stop_function {
    public:
      typedef T                                           value_type;
      typedef decltype(std::abs(T()))                     real_type;
      typedef std::complex<real_type>                     complex_type;
        
      typedef std::function<bool( CORK::vector<T>             const& eigvec, 
                                  real_type                   const& resid, 
                                  real_type                   const& norm_H, 
                                  T                           const& eigval)
                           >                              function_type;

    public:
      bool operator()( CORK::vector<value_type> const& eigen_vec, real_type const& resid, real_type norm_H, value_type const & val ) {
        return stop_function_(eigen_vec, resid, norm_H, val);
      }
    
    private:
      function_type     stop_function_;
      real_type         relative_tolerance_;
      real_type         absolute_tolerance_;
      
  };

  template <typename T>
  class stop_criterion_function {
    public:
      typedef T                                           value_type;
      typedef decltype(std::abs(T()))                     real_type;
      typedef std::complex<real_type>                     complex_type;

      typedef std::function<bool( CORK::vector<T>             const& eigvec, 
                                  real_type                   const& resid, 
                                  real_type                   const& norm_H, 
                                  T                           const& eigval)
                           >                              function_type;

    public:
      stop_criterion_function()
      : value_([]( CORK::vector<value_type> const& eigen_vec, real_type const& resid, real_type norm_H, value_type const & val) { return true; })
      {}

      template <typename ...Ts>
      stop_criterion_function( std::tuple<Ts...> const& options)
      : value_([]( CORK::vector<value_type> const& eigen_vec, real_type const& resid, real_type norm_H, value_type const & val) { return true; })
      {}

      // stop criterion based on the relative and absolute tolerance 
      stop_criterion_function( real_type relative_tol, real_type absolute_tol = 0. )
      : value_([=](CORK::vector<value_type> const& eigen_vec, real_type const& resid, real_type norm_H, value_type const & val) {
                      return resid < std::max(relative_tol*norm_H, absolute_tol);})
      {}

      // stop criterion based on the residual of the nep
      template <typename NEP>
      stop_criterion_function( NEP const& nep, real_type absolute_tol )
      : residu_vector(nep.size())
      , value_([&](CORK::vector<value_type> const& eigen_vec, real_type const& resid, real_type norm_H, value_type const & val) {
                      // beter om mee te geven via constructie van stop criterion zodat die niet altijd
                      // opnieuw gealloceerd wordt?
                      fill(residu_vector, value_type(0.));
                      nep.multiply_add( val, eigen_vec, residu_vector);
                      // kan hier relatieve tol gebruikt worden?
                      return norm_2(residu_vector) < absolute_tol;})
      {}

      // stop criterion given as a function by the user
      stop_criterion_function( function_type stopfunction )
      : value_(stopfunction)
      {}

      function_type  value() const { return value_ ; }
      function_type& value()       { return value_ ; }

    private:
      function_type               value_ ;
      CORK::vector<value_type>    residu_vector; // enkel indien nodig
      real_type                   relative_tolerance;
      real_type                   absolute_tolerance;
      bool                        needs_eigen_vec;
  } ;

} } // CORK::options

#endif
