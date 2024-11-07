//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_potential_theory_approximation_hpp
#define cork_approximation_potential_theory_approximation_hpp

#include <cork/basis/rational_newton.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cmath>
#include <limits>

namespace CORK { namespace approximation {

  template <typename T>
  class potential_theory_approximation {
  public:
    typedef T                       value_type ;
    typedef decltype(std::abs(T())) real_type ;

  private:
    typedef glas2::shared_vector<T> vector_type ;
    typedef glas2::shared_matrix<T> matrix_type ;

  public:
    typedef typename vector_type::size_type size_type ;

  public:
    potential_theory_approximation(size_type grade, size_type n_fun)
    : basis_( vector_type(grade+1), vector_type(grade), vector_type(grade+1) )
    , coefficients_(grade+1, n_fun)
    {}

  public:
    int size() const { return basis_.grade() ; }
    auto const& basis() const { return basis_ ; }
    matrix_type const& coefficients() const { return coefficients_ ; }
    real_type& error() { return error_ ;}
    real_type const& error() const { return error_ ;}

    matrix_type& coefficients() { return coefficients_ ; }
    vector_type& nodes() { return basis_.nodes() ; }
    vector_type& poles() { return basis_.poles() ; }
    vector_type& scaling() { return basis_.scaling() ; }

    vector_type const& nodes() const { return basis_.nodes() ; }
    vector_type const& poles() const { return basis_.poles() ; }
    vector_type const& scaling() const { return basis_.scaling() ; }

  public:
    template <typename Vector>
    void eval( value_type const& z, Vector result ) {
      assert( result.size()==coefficients_.num_columns() ) ;
      vector_type coefs( coefficients_.num_rows() ) ;
      basis_.evaluate( z, coefs ) ;
      result = multiply( transpose(coefficients_), coefs ) ;
    } // eval()

  private:
    basis::rational_newton< vector_type > basis_ ;
    matrix_type                           coefficients_ ;
    real_type                             error_ ;
  } ;

} } // namespace CORK::approximation

#endif
