//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_type_seed_hpp
#define glas2_type_seed_hpp

#include <random>
#include <chrono>
#include <type_traits>
#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 {

  template <typename T, typename EnableIf=void>
  class seed {
  } ;

  template <typename T>
  class seed< T, typename std::enable_if< std::is_floating_point<T>::value >::type > {
    public:
      seed()
      : distribution_( -1.0, 1.0 )
      {
        typedef std::chrono::system_clock myclock;
        myclock::duration d = myclock::now()- /*myclock::now() ; */myclock::time_point::min() ;
        unsigned seed2 = d.count();
        generator_.seed( seed2 ) ;
      }

    public:
      template <typename T2>
      typename std::enable_if< std::is_floating_point<T2>::value >::type var_gen( T2& t) const {
        t = T2( distribution_(generator_) ) ;
      }

#ifdef GLAS_COMPLEX
      template <typename T2>
      typename std::enable_if< std::is_floating_point<T2>::value >::type var_gen( std::complex<T2>& t) const {
        t.real( distribution_(generator_) ) ;
        t.imag( distribution_(generator_) ) ;
      }
#endif

    public:
      // We can use seed as a numerical type
      typedef seed<T> value_type ;

      operator T() const {
        return distribution_(generator_) ;
      }

#ifdef GLAS_COMPLEX
      operator std::complex<T>() {
        T temp_r, temp_i ;
        var_gen( temp_r ) ;
        var_gen( temp_i ) ;
        return std::complex<T>(temp_r,temp_i) ;
      }
#endif

      template <typename TT>
      operator TT() {
        TT val ;
        T temp ;
        var_gen( temp ) ;
        val = temp ;
      }


    private:
      mutable std::uniform_real_distribution<T> distribution_ ;
      mutable std::default_random_engine        generator_ ;
  } ;


} // namespace glas2

#endif
