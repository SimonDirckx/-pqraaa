//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_Q_DROP_TOL_HPP
#define CORK_OPTIONS_Q_DROP_TOL_HPP

#include<tuple>
#include<limits>

namespace CORK { namespace options {

  template <typename T>
  class q_drop_tol
  {
    public:
      q_drop_tol()
      : value_( std::numeric_limits<T>::epsilon() )
      {}

      template <typename ...Ts>
      q_drop_tol( std::tuple<Ts...> const& options )
      : value_( std::numeric_limits<T>::epsilon() )
      {}

      q_drop_tol( T const& v )
      : value_(v)
      {}

      T value() const { return value_ ; }
      T& value() { return value_ ; }

    private:
      T value_ ;
  } ;



} } // CORK::options

#endif
