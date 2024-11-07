//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_nonlinear_hpp
#define cork_matrix_valued_function_nonlinear_hpp

#include <cork/basis/union_of_functions.hpp>
#include <cork/matrix_valued_function/matrix_polynomial.hpp>

namespace CORK { namespace matrix_valued_function {

  template <typename Basis, typename FunctionSequence, typename CoefficientMatrices>
  decltype (auto) make_nonlinear( Basis const& basis, FunctionSequence const& function_sequence, CoefficientMatrices const& matrices )  {
    return make_matrix_polynomial( basis::union_of_functions<Basis,FunctionSequence>( basis, function_sequence ), matrices ) ;
  } // make_nonlinear()

  template <typename Basis, typename FunctionSequence, typename CoefficientMatrices>
  decltype (auto) make_nonlinear_lvalue( Basis const& basis, FunctionSequence const& function_sequence, CoefficientMatrices const& matrices )  {
    return make_matrix_polynomial( basis::union_of_functions<Basis const&,FunctionSequence const&>( basis, function_sequence ), matrices ) ;
  } // make_nonlinear_lvalue()

/*
  template <typename Basis, typename FunctionSequence, typename CoefficientMatrices, typename T=typename std::common_type< typename CoefficientMatrices::value_type, typename FunctionSequence::value_type>::type >
  class nonlinear
  {
    public:
      typedef typename std::decay<Basis>::type                basis_type ;
      typedef typename std::decay<FunctionSequence>::type     functions_sequence_type ;
      typedef typename std::decay<CoefficientMatrices>::type  matrix_type ;
      typedef T                                               value_type ;

    public:
      nonlinear( basis_type const& basis, functions_sequence_type const& functions_sequence, matrix_type const& matrices )
      : basis_( basis )
      , functions_sequence_( functions_sequence )
      , matrices_( matrices )
      {
        assert( matrices_.num_matrices() == basis_.num_terms()+functions_sequence_.num_terms() ) ;
      } // nonlinear

    public:
      auto num_terms() const { return basis_.num_terms() + functions_sequence_.num_terms() ; } // num_terms()
      auto size() const {
        assert( coefficient_matrices().num_rows() == coefficient_matrices().num_columns() ) ;
        return coefficient_matrices().num_rows() ;
      }

      basis_type const& basis() const { return basis_ ; }
      functions_sequence_type const& functions_sequence() const { return functions_sequence_; }
      matrix_type const& coefficient_matrices() const { return matrices_ ; }

    public:
      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w ) const {
        glas2::vector< typename basis_type::template value_type<Shift> > coefs( coefficient_matrices().num_matrices() ) ;
        basis_.evaluate( shift, coefs(glas2::range(0,basis_.num_terms())) ) ;
        functions_sequence_.evaluate( shift, coefs(glas2::range_from_end(basis_.num_terms(),0)) ) ;

        // Can lose efficiency. To be fixed.
        for (typename decltype(coefs)::size_type i=0; i<coefs.size(); ++i)
          matrices_.multiply_add( i, coefs(i)*x, w ) ;
      } // multiply_add()

    private:
      Basis                basis_ ;
      FunctionSequence     functions_sequence_ ;
      CoefficientMatrices  matrices_ ;
  } ; // class nonlinear

  template <typename Basis, typename FunctionSequence, typename CoefficientMatrices>
  decltype (auto) make_nonlinear_lvalue( Basis const& basis, FunctionSequence const& function_sequence, CoefficientMatrices const& matrices )  {
    return nonlinear< Basis const&, FunctionSequence const&, CoefficientMatrices const& >( basis, function_sequence, matrices ) ;
  }

  template <typename Basis, typename FunctionSequence, typename CoefficientMatrices>
  decltype (auto) make_nonlinear( Basis const& basis, FunctionSequence const& function_sequence, CoefficientMatrices const& matrices )  {
    return nonlinear< Basis, FunctionSequence, CoefficientMatrices >( basis, function_sequence, matrices ) ;
  }
*/
} } // namespace CORK::matrix_valued_function

#endif
