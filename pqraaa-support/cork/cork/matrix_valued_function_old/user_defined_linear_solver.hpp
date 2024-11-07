//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_user_defined_linear_solver_hpp
#define cork_matrix_valued_function_user_defined_linear_solver_hpp

#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/matrix_valued_function/linear_solver.hpp>
#include <cork/matrix_valued_function/nonlinear.hpp>
#include <cork/coefficient_matrices_linear_solver/user_defined.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace matrix_valued_function {

  template <typename Basis, typename Shift, typename CoefficientMatrices>
  class linear_solver< Basis, ::CORK::coefficient_matrices_linear_solver::linear_solver< Shift, coefficient_matrices::user_defined<CoefficientMatrices> > >
  {
    public:
      typedef coefficient_matrices_linear_solver::linear_solver< Shift, coefficient_matrices::user_defined<CoefficientMatrices> > linear_solver_type ;
      typedef typename std::decay< Basis >::type                           sequence_of_functions_type ;
      typedef typename linear_solver_type::size_type                       size_type ;

    public:
      typedef typename linear_solver_type::value_type value_type ;

    public:
      explicit linear_solver( Basis basis, linear_solver_type const& solver )
      : basis_( basis )
      , solver_( solver )
      {
        std::cout << "User defined" << std::endl ;
      }

    private:

    public:
      size_type num_terms() const { return basis_.num_terms() ; }

      size_type size() const { return solver_.size() ; }

   private:
      // Old codes
      void prepare_solve( value_type const& shift ) {
        // This is different from the way it is done for Dense and Sparse in CORK.
        solver_.prepare_solve( shift ) ;
      } // prepare_solve()

      template <typename V>
      void solve( V v ) const {
        static_assert( std::is_same< typename std::decay<V>::type::value_type, value_type >::value, "V must have the same value_type as solver" ) ;
        solver_.solve( v ) ;
      } // solve()

    public:
      // For 2020 code
      template <typename V, typename Options>
      void solve( value_type& shift, V v, bool is_new_shift, Options const& options ) const {
        solver_.solve( shift, v, is_new_shift ) ;
      }

    public:
      typedef linear_solver< sequence_of_functions_type const&, typename linear_solver_type::clone_type > clone_type ;
      clone_type clone() const {
        return clone_type( basis_, solver_.clone() ) ;
      }

    private:
      Basis              basis_ ;
      linear_solver_type solver_ ;
  } ; // linear_solver


  template <typename ValueType, typename MatrixPolynomial, typename CoefficientMatrices>
  struct linear_solver_struct< ValueType, MatrixPolynomial, coefficient_matrices::user_defined<CoefficientMatrices> > {
    typedef typename MatrixPolynomial::template value_type<ValueType>                                                                       value_type ;
    typedef ::CORK::coefficient_matrices_linear_solver::linear_solver< ValueType, coefficient_matrices::user_defined<CoefficientMatrices> > solver_type ;
    typedef linear_solver< typename MatrixPolynomial::basis_type const&, solver_type >                                                      type ;

    static type apply( MatrixPolynomial const& mp ) {
      std::cout << "Select this one" << std::endl ;
      return type( mp.basis(), solver_type( mp.coefficient_matrices() ) ) ;
    }
  } ;


} } // namespace CORK::matrix_valued_function

#endif
