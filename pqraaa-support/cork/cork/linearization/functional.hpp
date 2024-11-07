//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_functional_hpp
#define cork_linearization_functional_hpp

#include <glas2/matrix.hpp>
#include <cassert>
#include <type_traits>

namespace CORK {

  template <typename PolynomialBasis, typename Matrix>
  class functional
  : public PolynomialBasis
  {
    public:
      typedef typename PolynomialBasis::value_type       value_type ;
      typedef typename PolynomialBasis::shift_value_type shift_value_type ;

      typedef PolynomialBasis                       polynomial_basis_type ;
      typedef Matrix                                functional_matrix_type ;
      typedef std::vector< functional_matrix_type > functional_matrix_sequence_type ;

    public:
      functional( polynomial_basis_type poly, functional_matrix_sequence_type const& matrices )
      : polynomial_basis_type( poly )
      , matrices_( matrices )
      , factor_( matrices[0].num_rows(), matrices[0].num_rows() )
      , pivot_( matrices_[0].num_rows() )
      {}

    public:
      void shift( shift_value_type s ) {
        static_cast<polynomial_basis_type&>(*this).shift(s) ;
        factor_ = ...
        int info = boost::numeric::bindings::lapack::getrs( factor_, pivot_, w ) ;
        assert( info>=0 ) ;
        if (info>0) throw std::runtime_error( "Matrix is singular" ) ;
      }

      template <typename QQ, typename UU, typename ZZ, typename W>
      void solve_P( QQ const& Q, UU const& U, ZZ const& z, W w ) const {
        w = ...
        boost::numeric::bindings::lapack::getrs( factor_, pivot_, w ) ;
        assert( info!=0 ) ;
      } // solve_P()

    private:
      functional_matrix_type matrices_ ;
      functional_matrix_type factor_ ;
      std::vector<int>  pivots_ ;
  } ;

} // namespace CORK

#endif
