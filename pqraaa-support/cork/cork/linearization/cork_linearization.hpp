//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_cork_linearization_hpp
#define cork_linearization_cork_linearization_hpp

#include <cork/utility/value_type_for.hpp>
#include <cork/linearization/cork_linearization_multiply_handle.hpp>
#include <cork/linearization/cork_linearization_solve_handle.hpp>
#include <cork/linearization/cork_linearization_fill_handle.hpp>
#include <cork/linearization/info.hpp>
#include <cork/matrix_valued_function/nonlinear_matrix.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/matrix_iterator/matrix_iterator.hpp>
#include <cork/utility/ref.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <type_traits>

namespace CORK { namespace linearization {

  template <typename MatrixPolynomial>
  class cork_linearization
  {
    public:
      typedef typename CORK::deref_type< MatrixPolynomial >::type        matrix_polynomial_type ;
      typedef typename matrix_polynomial_type::basis_type                basis_type ;
      typedef typename matrix_polynomial_type::coefficient_matrices_type coefficient_matrices_type ;
      typedef basis4CORK::basis4CORK< basis_type >                       basis4cork_type ;
      typedef matrix_iterator::matrix_iterator< basis_type >             matrix_iterator_type ;

      typedef typename basis4cork_type::size_type                        size_type ;

    public:
      cork_linearization( MatrixPolynomial poly, info& information )
      : matrix_polynomial_( poly )
      , basis_( CORK::deref(matrix_polynomial_).basis() )
      , matrix_iterator_( CORK::deref(matrix_polynomial_).basis() )
      , information_( information )
      {
        information_.linearization = "CORK Linearization" ;
      }

    public:
      size_type num_rows() const { return size_of_basis()*size() ; }
      size_type num_columns() const { return size_of_basis()*size() ; }

      size_type size_of_basis() const { return basis_.size() ; }
      size_type size() const {
        return CORK::deref(matrix_polynomial_).size() ;
      }

      auto const& matrix_polynomial() const { return CORK::deref(matrix_polynomial_) ; }

      basis4cork_type const& basis() const { return basis_ ; }
      coefficient_matrices_type const& coefficient_matrices() const { return CORK::deref(matrix_polynomial_).coefficient_matrices() ; }
      matrix_iterator_type const& matrix_iterator() const { return matrix_iterator_ ; }
      matrix_iterator_type& matrix_iterator() { return matrix_iterator_ ; }
      info& information() { return information_ ; }

    /*public:
      typedef cork_linearization< basis4cork_type, matrix_iterator_type&, coefficient_matrices_type const&, typename linear_solver_type::clone_type > clone_type ;
      clone_type clone() {
        return clone_type( basis_, matrix_iterator_, coefficient_matrices_, linear_solver_.clone(), information_ ) ;
      }*/

  public:
    template <typename T>
    using value_type_for = typename matrix_polynomial_type::template value_type_for<T> ;

    public:
      template <typename T>
      using multiply_handle_type = cork_linearization_multiply_handle< T
                                                                     , basis4cork_type const&
                                                                     , matrix_iterator_type const&
                                                                     , coefficient_matrices_type const&
                                                                     > ;

      template <typename T>
      decltype(auto) multiply_handle() const {
        return multiply_handle_type<T>( basis_, matrix_iterator_, coefficient_matrices(), information_ ) ;
      }

    public:
      using fill_handle_type = cork_linearization_fill_handle< basis4cork_type const&
                                                             , matrix_iterator_type const&
                                                             , coefficient_matrices_type const&
                                                             > ;

      decltype(auto) fill_handle() const {
        return fill_handle_type( basis_, matrix_iterator_, coefficient_matrices() ) ;
      }

    public:
      template <typename Solver>
      using solve_handle_external_type = cork_linearization_solve_handle< Solver
                                                                        , basis4cork_type const&
                                                                        , typename basis4cork_type::template handle_type<typename std::decay<Solver>::type::value_type>
                                                                        , matrix_iterator_type const&
                                                                        , coefficient_matrices_type const&
                                                                        > ;

      template <typename Solver>
      decltype(auto) solve_handle_external_external( Solver& linear_solver) const {
        typedef typename Solver::value_type solver_value_type ;
        return solve_handle_external_type<Solver&>( basis_
                                                  , basis_.template handle<solver_value_type>()
                                                  , matrix_iterator_
                                                  , coefficient_matrices()
                                                  , linear_solver
                                                  , information_
                                                  ) ;
      }

      template <typename ValueType>
      using solve_handle_type = cork_linearization_solve_handle< typename matrix_polynomial_type::template linear_solver_type<ValueType>
                                                               , basis4cork_type const&
                                                               , typename basis4cork_type::template handle_type<typename matrix_polynomial_type::template linear_solver_type<ValueType>::value_type>
                                                               , matrix_iterator_type const&
                                                               , coefficient_matrices_type const&
                                                               > ;


      template <typename ShiftValueType>
      decltype(auto) solve_handle() const {
        auto linear_solver = CORK::deref(matrix_polynomial_).template linear_solver< ShiftValueType >() ; //matrix_valued_function::make_linear_solver<ShiftValueType>( matrix_polynomial_ ) ;
        typedef typename decltype(linear_solver)::value_type solver_value_type ;
        return solve_handle_type<ShiftValueType>( basis_
                                                , basis_.template handle<solver_value_type>()
                                                , matrix_iterator_
                                                , coefficient_matrices()
                                                , linear_solver
                                                , information_ 
                                                ) ;
      }

    private:
      MatrixPolynomial                 matrix_polynomial_ ;
      basis4cork_type                  basis_ ;
      matrix_iterator_type             matrix_iterator_ ;
      info&                            information_ ;
  } ; // default_linearization


  template <typename MatrixPolynomial>
  auto make_cork_linearization( MatrixPolynomial const& mp, info& information ) {
    return cork_linearization< MatrixPolynomial >( mp, information ) ;
  }

} } // namespace CORK::linearization


namespace CORK {

  template <typename T, typename MatrixPolynomial>
  struct value_type_for< T, linearization::cork_linearization< MatrixPolynomial > >
  {
    typedef typename MatrixPolynomial::template value_type_for<T> type ;
  } ;

} // namespace CORK

#endif
