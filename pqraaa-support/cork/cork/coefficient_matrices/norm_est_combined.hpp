//  (C) Copyright Karl Meerbergen 2022.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompany_glasing file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_norm_est_combined_hpp
#define cork_coefficient_matrices_norm_est_combined_hpp

#include <cork/coefficient_matrices/norm_est.hpp>

namespace CORK { namespace coefficient_matrices {

  template <typename T>
  class norm_est_combined
  : public norm_est< T >
  {
    private:
      typedef norm_est< T > base_type ;

    public:
      template <typename NormEst, typename Combined>
      norm_est_combined( NormEst const& norm_est, Combined const& C )
      : base_type( norm_est.R().num_rows(), C.num_columns() )
      //, n_( C.num_rows() )
      {
        assert( C.num_rows()==norm_est.R().num_columns() ) ;

        this->R_ = multiply( norm_est.R(), C ) ;
      }
  } ; // norm_est_combined

} } // namespace CORK::coefficient_matrices

#endif
