//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_nep_residual_norm_hpp
#define cork_eigs_nep_residual_norm_hpp

#include <glas2/vector.hpp>

namespace CORK { namespace eigs {

 template <typename MatrixValuedFunction, typename EigenValue, typename EigenVector, typename Temp, typename Resid>
 void nep_residual_norm( MatrixValuedFunction const& nep, EigenValue const& eigen_value, EigenVector const& eigen_vector, Temp temp_vector1, Resid& resid ) {
   assert( nep.size()==eigen_vector.size() ) ;
   assert( temp_vector1.size()==eigen_vector.size() ) ;

   fill(temp_vector1,0.0) ;
   nep.multiply_add( eigen_value, eigen_vector, temp_vector1 ) ;
   resid = norm_2( temp_vector1 ) ;
 } // nep_residual_norm()
   
} } // namespace CORK::eigs


#endif
