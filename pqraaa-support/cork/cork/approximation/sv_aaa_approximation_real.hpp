//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_approximation_real_hpp
#define cork_approximation_sv_aaa_approximation_real_hpp

#include <cork/approximation/sv_aaa_approximation.hpp>

namespace CORK { namespace approximation {

  template <typename T>
  class SV_AAA_approximation_real
  : public SV_AAA_approximation<T,T,decltype(std::abs(T()))>
  {
  private:
    typedef decltype(std::abs(T())) real_type ;
    typedef SV_AAA_approximation<T,T,real_type> base_type ;

  public:
    SV_AAA_approximation_real(typename base_type::size_type n_fun, typename base_type::size_type n_buf = 100)
    : base_type( n_fun, n_buf )
    {}

    template <typename Vector>
    void eval( T const& z, Vector result ) const {
      auto nodes = this->nodes() ;
      auto weights = this->weights() ;
      auto coefficients = this->coefficients() ;
      glas2::vector< T > temp( result.size() ) ;
      T sum = 0.0 ;
      fill( result, 0.0 ) ;
      for (typename base_type::size_type i=0; i<this->n(); ++i) {
        if (nodes(i).imag()==0.0) {
          if (z==nodes(i)) {
            result = coefficients(i, glas2::all()) ;
            return ;
          }
          if (nodes(i)==std::numeric_limits<real_type>::infinity()) {
            result += weights(i) * coefficients(i, glas2::all()) ;
            sum += weights(i) ;
          } else {
            result += weights(i) * coefficients(i, glas2::all()) / (z - nodes(i) ) ;
            sum += weights(i) / (z - nodes(i) ) ;
          }
        } else {
          assert( nodes(i)!=std::numeric_limits<real_type>::infinity() ) ;
          if (z==nodes(i)) {
            result = coefficients(i, glas2::all()) ;
            result += T(0.,1.0)*coefficients(i+1, glas2::all()) ;
            return ;
          } else if (z==conj(nodes(i))) {
            result = coefficients(i, glas2::all()) ;
            result -= T(0.,1.0)*coefficients(i+1, glas2::all()) ;
            return ;
          }
          temp = coefficients(i, glas2::all()) ;
          temp += T(0.,1.0)*coefficients(i+1, glas2::all()) ;

          result += weights(i) * temp / (z - nodes(i) ) ;
          result += weights(i+1) * conj(temp) / (z - nodes(i+1) ) ;
          sum += weights(i) / (z - nodes(i) ) + weights(i+1) / (z - nodes(i+1) ) ;
          ++i ;
        }
      }
      if (this->n()>0) {
   /*     if (sum==0.0) {
          // Infinite or very large z. Multiply all terms by 1/z
          for (typename base_type::size_type i=0; i<this->n(); ++i) {
            if (nodes(i).imag()==0.0) {
              assert( z!=nodes(i)) ;
              if (nodes(i)==std::numeric_limits<real_type>::infinity()) {
                result += weights(i) * coefficients(i, glas2::all()) * z ;
                sum += weights(i) * z ;
              } else {
                result += weights(i) * coefficients(i, glas2::all()) / (z - nodes(i) ) ;
                sum += weights(i) / (z - nodes(i) ) ;
              }
            } else {
              assert( nodes(i)!=std::numeric_limits<real_type>::infinity() ) ;
              assert( z!=nodes(i)) ;
              assert( z!=conj(nodes(i))) ;
              temp = coefficients(i, glas2::all()) ;
              temp += T(0.,1.0)*coefficients(i+1, glas2::all()) ;

              result += weights(i) * temp / (z - nodes(i) ) ;
              result += weights(i+1) * conj(temp) / (z - nodes(i+1) ) ;
              sum += weights(i) / (z - nodes(i) ) + weights(i+1) / (z - nodes(i+1) ) ;
              ++i ;
            }
          }
        }*/
        result /= sum ;
      }
    } // eval()
  } ; // class SV_AAA_approximation

} } // namespace CORK::approximation

#endif
