//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_AAA_STOP_CRITERION_HPP
#define CORK_OPTIONS_AAA_STOP_CRITERION_HPP

#include <limits>
#include <tuple>

namespace CORK { namespace options {

  template <typename T>
  class aaa_stop_criterion {
    public:
      enum enum_type { MAX, WEIGHTED } ;

    public:
      aaa_stop_criterion()
      : value_( enum_type::MAX )
      , tolerance_( std::numeric_limits<T>::epsilon() )
      {}

      template <typename ...Ts>
      aaa_stop_criterion( std::tuple<Ts...> const& options)
      : value_( enum_type::MAX )
      , tolerance_( std::numeric_limits<T>::epsilon() )
      {}

      aaa_stop_criterion( enum_type const& v, T const& tolerance )
      : value_(v)
      , tolerance_( tolerance )
      {}

      auto const& value() const { return *this ; }

      T const& tolerance() const { return tolerance_ ; }
      enum_type choice() const { return value_ ;}

    public:
      static aaa_stop_criterion max() { return aaa_stop_criterion( MAX ) ; }
      static aaa_stop_criterion weighted() { return aaa_stop_criterion( WEIGHTED ) ; }

    private:
      enum_type     value_ ;
      T             tolerance_ ;
  } ;


  template <typename T>
  auto aaa_stop_criterion_max( T const& tolerance=0. ) {
    return aaa_stop_criterion<T>( aaa_stop_criterion<T>::MAX, tolerance ) ;
  }

  template <typename T>
  auto aaa_stop_criterion_weighted( T const& tolerance=0. ) {
    return aaa_stop_criterion<T>( aaa_stop_criterion<T>::WEIGHTED, tolerance ) ;
  }

} } // CORK::options

#endif
