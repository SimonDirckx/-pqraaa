//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_eigen_residual_norms_hpp
#define cork_eigs_eigen_residual_norms_hpp

#include <cork/eigs/eigen_pairs.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace eigs {

 template <typename MatrixValuedFunction, typename EigenPairs, typename Resid>
 void eigen_residual_norms( MatrixValuedFunction const& nep, EigenPairs const& eig_pairs, Resid resid ) {
   typedef typename EigenPairs::value_type value_type ;
   assert( nep.size()==eig_pairs.vectors().num_rows() ) ;
   assert( eig_pairs.size()==resid.size() ) ;

   glas2::shared_vector< value_type > res( nep.size() ) ;

   for (int i=0; i<eig_pairs.size(); ++i) {
     fill(res,0.0) ;
     nep.multiply_add( eig_pairs.eigenvalues()(i), eig_pairs.vectors()(glas2::all(),i), res ) ;
     resid(i) = norm_2( res ) ;
   }
 } // eigen_residual_norms()


 template <typename MatrixValuedFunction, typename EigenValues, typename EigenVectors, typename Resid>
 void eigen_residual_norms( MatrixValuedFunction const& nep, EigenValues const& eig_values, EigenVectors const& eig_vectors, Resid resid ) {
   typedef typename EigenValues::value_type value_type ;
   assert( nep.size()==eig_vectors.num_rows() ) ;
   assert( eig_values.size()==resid.size() ) ;

   glas2::shared_vector< value_type > res( nep.size() ) ;

   for (int i=0; i<eig_values.size(); ++i) {
     fill(res,0.0) ;
     nep.multiply_add( eig_values(i), eig_vectors(glas2::all(),i), res ) ;
     resid(i) = norm_2( res ) ;
   }
 } // eigen_residual_norms()
   
} } // namespace CORK::eigs


#endif
