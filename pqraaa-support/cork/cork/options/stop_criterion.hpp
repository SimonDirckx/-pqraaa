//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)
#ifndef CORK_OPTIONS_STOP_CRITERION_HPP
#define CORK_OPTIONS_STOP_CRITERION_HPP

#include <cork/vector.hpp>
#include <functional>
#include <limits>
#include <tuple>

namespace CORK { namespace options {

  template <typename T>
  class stop_criterion {
    public:
      typedef T                                           value_type;
      typedef decltype(std::abs(T()))                     real_type;
      typedef std::complex<real_type>                     complex_type;

      enum enum_type { KRYLOV, NEP, MIXED, INVARIANT_PAIR } ;

    public:
      stop_criterion()
      : choice_( KRYLOV )
      , relative_tolerance_(0.0)
      , absolute_tolerance_(0.0)
      {}

      stop_criterion( enum_type choice, real_type const& relative_tolerance, real_type const& absolute_tolerance )
      : choice_( choice )
      , relative_tolerance_(relative_tolerance)
      , absolute_tolerance_(absolute_tolerance)
      {}

      template <typename ...Ts>
      stop_criterion( std::tuple<Ts...> const& options)
      : choice_( KRYLOV )
      , relative_tolerance_(0.0)
      , absolute_tolerance_(0.0)
      {}

    public:
      // stop criterion based on the relative and absolute tolerance 
      auto value() const { return *this ; }
      auto& value()       { return *this ; }

    public:
      real_type const& relative_tolerance() const { return relative_tolerance_ ; }
      real_type const& absolute_tolerance() const { return absolute_tolerance_ ; }
      enum_type choice() const { return choice_ ; }

    private:
      enum_type       choice_ ;
      real_type       relative_tolerance_;
      real_type       absolute_tolerance_;
  } ; // class stop_criterion


  template <typename RealType>
  auto stop_criterion_rks( RealType relative_tol=0.0, RealType absolute_tol = 0. ) {
    return stop_criterion<RealType>( stop_criterion<RealType>::KRYLOV, relative_tol, absolute_tol ) ;
  }

  template <typename RealType>
  auto stop_criterion_nep( RealType relative_tol=0.0, RealType absolute_tol = 0. ) {
    return stop_criterion<RealType>( stop_criterion<RealType>::NEP, relative_tol, absolute_tol ) ;
  }

  template <typename RealType>
  auto stop_criterion_mixed( RealType relative_tol=0.0, RealType absolute_tol = 0. ) {
    return stop_criterion<RealType>( stop_criterion<RealType>::MIXED, relative_tol, absolute_tol ) ;
  }

  template <typename RealType>
  auto stop_criterion_invariant_pair( RealType relative_tol=0.0, RealType absolute_tol = 0. ) {
    return stop_criterion<RealType>( stop_criterion<RealType>::INVARIANT_PAIR, relative_tol, absolute_tol ) ;
  }

} } // CORK::options

#endif
