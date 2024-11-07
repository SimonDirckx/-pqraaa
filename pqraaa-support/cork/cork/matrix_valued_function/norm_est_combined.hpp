//  (C) Copyright Karl Meerbergen 2022.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_norm_est_combined_hpp
#define cork_matrix_valued_function_norm_est_combined_hpp

#include <cork/coefficient_matrices/norm_est_combined.hpp>
#include <cork/matrix_valued_function/norm_est.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {


  template <typename NormEst, typename Combination>
  class norm_est_combined
  : public deref_type<NormEst>::type
  {
    public:
      typedef typename deref_type<NormEst>::type  base_type ;

    private:
      typedef coefficient_matrices::norm_est_combined< typename std::common_type< typename base_type::value_type, typename std::decay<Combination>::type::value_type >::type > norm_est_type ;
      //typedef coefficient_matrices::norm_est_combined< typename base_type::value_type > norm_est_type ;

    public:
      norm_est_combined( base_type const& norm_est, Combination combination  )
      : base_type( norm_est )
      , combination_( combination )
      {}

    public:
      void activate() {
        static_cast<base_type&>(*this).activate() ;
#ifndef NDEBUG
        if (norm_est_combined_.get()) std::cerr << "WARNING: CORK::matrix_valued_function::norm_est_combined::activate(): is already activated" << std::endl ;
#endif
        norm_est_combined_.reset( new norm_est_type( this->coefficients_norm_est(), combination_ ) ) ;
        //norm_est_combined_.reset( new norm_est_type( *this, combination_ ) ) ;
      }

    public:
      auto norm_coefficient( int i ) const {
        assert( norm_est_combined_.get() ) ;
        return norm_est_combined_->norm_coefficient(i) ;
      } // norm_coefficient() ;

    private:
      std::shared_ptr<norm_est_type> norm_est_combined_ ;
      Combination                    combination_ ;
  } ; // class norm_est_combined

} } // namespace CORK::matrix_valued_function

#endif
