//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_set_valued_function_hpp
#define cork_approximation_set_valued_function_hpp

#include <glas2/vector.hpp>
#include <cmath>
#include <limits>

namespace CORK { namespace approximation {

  template <typename T, typename Basis>
  class set_valued_function {
  public:
    typedef T                       value_type ;

  public:
    typedef glas2::shared_matrix<value_type> matrix_type ;

  public:
    typedef typename Basis::size_type size_type ;

  public:
    template <typename C>
    set_valued_function( Basis const& basis, C const& coefficients )
    : basis_( basis )
    , coefficients_(coefficients.num_rows(), coefficients.num_columns())
    {
      assert( basis.num_terms()==coefficients.num_rows() ) ;
      coefficients_ = coefficients ;
    }

  public:
    size_type n() const { return basis_.num_terms() ; }

  public:
    auto const& basis() const { return basis_ ; }
    auto const& coefficients() const { return coefficients_ ; }
    auto& coefficients() { return coefficients_ ; }

  public:
    template <typename Arg, typename Vector>
    void eval( Arg const& z, Vector result ) const {
      assert( result.size() == coefficients_.num_columns() ) ;

      glas2::vector< typename Basis::template value_type_for<Arg> > coefs( basis_.num_terms() ) ;
      basis_.evaluate( z, coefs ) ;

      result = multiply( transpose(coefficients_), coefs ) ;
    } // eval()

  private:
    Basis       basis_ ;
    matrix_type coefficients_ ;
  } ;

} } // namespace CORK::approximation

#endif
