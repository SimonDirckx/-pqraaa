//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_BACKEND_HPP
#define CORK_OPTIONS_BACKEND_HPP

#include <cork/backend/default_backend.hpp>
#include <cassert>
#include<tuple>

namespace CORK { namespace options {

  struct backend_key
  {
    public:
      backend_key()
      {}

      template <typename Options>
      backend_key( Options const& )
      {}

      auto value() const { return backend::default_backend() ; }
  } ;

  template <typename B>
  class backend_selector
  : public backend_key
  {
    public:
      backend_selector( B const& b )
      : backend_( &b )
      {}

      template <typename Options>
      backend_selector( Options const& )
      : backend_( 0 )
      {
        assert( false ) ;
      }

    public:
      B const& value() const { return *backend_ ; }

    private:
      B const* backend_ ;
  } ;

  template <typename B>
  backend_selector<B> backend( B const& b ) { return backend_selector<B>(b) ; }

} } // CORK::options

#endif
