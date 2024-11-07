//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_type_real_proxy_hpp
#define glas2_type_real_proxy_hpp

#include <glas2/concept/is_arithmetic.hpp>
#include <type_traits>
#include <cassert>
#ifdef GLAS_COMPLEX
#include <complex>
#endif

#ifdef glas2_concept_abs_squared_hpp
#error <glas2/type/real_proxy.hpp> should be included before <glas2/concept/abs_squared.hpp>
#endif

namespace glas2 {

  template <typename T>
  class real_proxy {
  public:
    typedef typename T::value_type value_type ;

    real_proxy()
#ifndef NDEBUG
    : value_(0)
#endif
    {}

    real_proxy( T& value )
    : value_( &value )
    {}

    operator value_type() const {
      assert( value_!=0 ) ;
      return value_->real() ;
    }

    template <typename E>
    real_proxy& operator=( E const& e ) {
      assert( value_!=0 ) ;
      value_->real( e ) ;
      //value_->real() = e ;
      return *this ;
    }

  private:
    T* value_ ;
  } ;

  template <typename T>
  struct is_arithmetic< real_proxy<T> >
  : is_arithmetic<T>
  {} ;

  template <typename T>
  typename std::enable_if< glas2::is_arithmetic<T>::value, typename T::value_type >::type abs_squared( real_proxy<T> const& x ) {
    return typename T::value_type(x) * typename T::value_type(x) ;
  }

} // namespace glas2

#endif
