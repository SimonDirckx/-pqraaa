//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_accumulated_partial_fractions_real_matrices_hpp
#define cork_basis4cork_accumulated_partial_fractions_real_matrices_hpp

#include <cork/utility/is_infinite.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>


namespace CORK { namespace basis4CORK {

  template <typename Weights, typename Poles, typename I, typename Matrix>
  void accumulated_partial_fractions_real_matrices( Weights const& weights, Poles const& poles, I num_real, Matrix M, Matrix N ) {
     fill(M, 0.0) ;
     fill(N, 0.0) ;

     if (num_real>0) {
       M(0,0) = weights(0).real() ;
     } else {
       M(0,0) = 2.*weights(0).real() ;
       M(1,0) = -2.*weights(0).imag() ;
     }
     typename Weights::value_type previous_weight = 1.0 ;

     // Real poles
     for (int i=1; i<=num_real; ++i) {
       // Diagonal block
       M(i-1,i) = previous_weight.real() * poles(i-1).real() ;
       N(i-1,i) = previous_weight.real() ;

       // Lower diagonal block
       if (i<M.num_rows()) {
         M(i,i) = weights(i).real() ;
       }
       previous_weight = weights(i-1) ;
     }

     // Transition real to complex poles
     if (num_real>0 && num_real<poles.size()) {
       // Diagonal block
       M(num_real-1,num_real) = previous_weight.real() * poles(num_real-1).real() ;
       N(num_real-1,num_real) = previous_weight.real() ;

       // Lower diagonal block
       M(num_real,num_real) = 2.0*weights(num_real).real() ;
       M(num_real+1,num_real) = -2.0*weights(num_real).imag() ;

       previous_weight = weights(num_real-1) ;
     }

     // Complex poles
     for (int i=num_real; i<M.num_rows(); i+=2) {
       // Diagonal block
       M(i,i+1) = std::real( poles(i)*previous_weight ) ;
       N(i,i+1) = previous_weight.real() ;

       M(i,i+2) = std::imag( poles(i)*previous_weight ) ;
       N(i,i+2) = previous_weight.imag() ;

       M(i+1,i+1) = -M(i,i+2) ; N(i+1,i+1) = -N(i,i+2) ;
       M(i+1,i+2) = M(i,i+1) ; N(i+1,i+2) = N(i,i+1) ;

       // Lower diagonal block
       if (i<M.num_rows()-2) {
         M(i+2,i+1) = std::real( weights(i+2) ) ;

         M(i+2,i+2) = std::imag( weights(i+2) ) ;

         M(i+3,i+1) = -M(i+2,i+2) ;
         M(i+3,i+2) = M(i+2,i+1) ;
       }

       previous_weight = weights(i) ;
     }
  } // accumulated_partial_fractions_real_matrices()

} } // namespace CORK::basis4cork

#endif
