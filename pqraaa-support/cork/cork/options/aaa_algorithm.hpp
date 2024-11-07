//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_AAA_ALGORITHM_HPP
#define CORK_OPTIONS_AAA_ALGORITHM_HPP

#include <limits>
#include <tuple>

namespace CORK { namespace options {

  class aaa_algorithm {
    public:
      enum enum_type { SV_AAA, SVD_SV_AAA } ;

    public:
      aaa_algorithm()
      : value_( enum_type::SVD_SV_AAA )
      {}

      template <typename ...Ts>
      aaa_algorithm( std::tuple<Ts...> const& options)
      : value_( enum_type::SVD_SV_AAA )
      {}

      aaa_algorithm( enum_type const& v )
      : value_(v)
      {}

      enum_type const&  value() const { return value_ ; }
      enum_type& value()       { return value_ ; }

    public:
      static aaa_algorithm sv_aaa() { return aaa_algorithm( SV_AAA ) ; }
      static aaa_algorithm svd_sv_aaa() { return aaa_algorithm( SVD_SV_AAA ) ; }

    private:
      enum_type     value_ ;
  } ;

} } // CORK::options

#endif
