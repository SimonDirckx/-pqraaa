//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_barycentric_rational_real_matrices_hpp
#define cork_basis4cork_barycentric_rational_real_matrices_hpp

#include <cork/utility/is_infinite.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>


namespace CORK { namespace basis4CORK {

  template <typename Weights, typename Nodes, typename Matrix>
  void barycentric_rational_real_matrices( Weights const& weights, Nodes const& nodes, Matrix M, Matrix N ) {
     fill(M, 0.0) ; M(0,0) = -1.0 ;
     fill(N, 0.0) ;

     int i_row = 1;
     for (int i=0; i_row<M.num_rows(); ++i) {
       if (nodes(i).imag()==0.0 && nodes(i+1).imag()==0.0) {
         M(0,i+1) = 1.0 ;
         if (is_infinite(nodes(i))) {
           M(i_row,i+1) = weights(i+1).real() ;

           M(i_row,i+2) = weights(i).real()*nodes(i+1).real() ;
           N(i_row,i+2) = weights(i).real() ;
           assert(!is_infinite(nodes(i+1))) ;
         } else {
           M(i_row,i+1) = -weights(i+1).real()*nodes(i).real() ;
           N(i_row,i+1) = -weights(i+1).real() ;

           if (is_infinite(nodes(i+1))) {
             M(i_row,i+2) = -weights(i).real() ;
           } else {
             M(i_row,i+2) = weights(i).real()*nodes(i+1).real() ;
             N(i_row,i+2) = weights(i).real() ;
           }
         }
         ++i_row ;
       } else if (nodes(i).imag()==0.0 && nodes(i+1).imag()!=0.0) {
         assert( !is_infinite(nodes(i+1)) ) ;
         M(0,i+1) = 1.0 ;
         if (is_infinite(nodes(i))) {
           M(i_row,i+1) = 2.0*weights(i+1).real() ;
           M(i_row+1,i+1) = -2.0*weights(i+1).imag() ;
         } else {
           M(i_row,i+1) = -2.0*weights(i+1).real()*nodes(i).real() ;
           N(i_row,i+1) = -2.0*weights(i+1).real() ;
           M(i_row+1,i+1) = 2.0*weights(i+1).imag()*nodes(i).real() ;
           N(i_row+1,i+1) = 2.0*weights(i+1).imag() ;
         }

         M(i_row,i+2) = weights(i).real()*nodes(i+1).real() ;
         N(i_row,i+2) = weights(i).real() ;
         M(i_row+1,i+2) = -weights(i).real()*nodes(i+1).imag() ;
         M(i_row,i+3) = -M(i_row+1,i+2) ;
         M(i_row+1,i+3) = M(i_row,i+2) ;
         N(i_row+1,i+3) = N(i_row,i+2) ;

         i_row += 2 ;
       } else /*if (nodes(i).imag()!=0.0)*/ {
         assert( nodes(i).imag()!=0.0 ) ;
         assert( !is_infinite(nodes(i)) ) ;
         if (i_row>=M.num_rows()-1) {
           M(0,i+1) = 1.0 ; M(0,i+2) = 0.0 ;

           auto w_n = nodes(i+1)*weights(i) ;
           M(i_row,i+1) = imag(w_n) ;
           N(i_row,i+1) = weights(i).imag() ;
           M(i_row,i+2) = real(w_n) ;
           N(i_row,i+2) = weights(i).real() ;
           ++i_row ;
         } else {
           assert( nodes(i+2).imag()!=0.0 ) ;

           M(0,i+1) = 1.0 ; M(0,i+2) = 0.0 ;

           auto w_n = nodes(i)*weights(i+2) ;
           M(i_row,i+1) = - real(w_n) ;
           N(i_row,i+1) = - weights(i+2).real() ;
           M(i_row,i+2) = -imag(w_n) ;
           N(i_row,i+2) = -weights(i+2).imag() ;
           M(i_row+1,i+1) = -M(i_row,i+2) ;
           N(i_row+1,i+1) = -N(i_row,i+2) ;
           M(i_row+1,i+2) = M(i_row,i+1) ;
           N(i_row+1,i+2) = N(i_row,i+1) ;

           w_n = nodes(i+2)*weights(i) ;
           M(i_row,i+3) = real(w_n) ;
           N(i_row,i+3) = weights(i).real() ;
           M(i_row,i+4) = imag(w_n) ;
           N(i_row,i+4) = weights(i).imag() ;
           M(i_row+1,i+3) = -M(i_row,i+4) ;
           N(i_row+1,i+3) = -N(i_row,i+4) ;
           M(i_row+1,i+4) = M(i_row,i+3) ;
           N(i_row+1,i+4) = N(i_row,i+3) ;

           ++i;
           i_row += 2 ;
         }
       }
     }
     if (nodes(nodes.size()-1).imag()==0.0) {
       M(0,nodes.size()) = 1.0 ;
     } else {
       M(0,nodes.size()-1) = 1.0 ;
       M(0,nodes.size()) = 0.0 ;
     }
  } // barycentric_rational_real_matrices()

} } // namespace CORK::basis4cork

#endif
