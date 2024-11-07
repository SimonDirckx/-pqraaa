//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_empty_hpp
#define cork_matrix_valued_function_empty_hpp

#include <cork/vector.hpp>
#include <cork/coefficient_matrices/matrices_by_functions.hpp>
#include <cork/matrix_valued_function/matrix_polynomial.hpp>
#include <cork/basis/set_of_functions.hpp>
#include <cork/exception/not_implemented.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {

  template <typename Arg>
  class empty {
    public:
      struct multiply_add_type {
        template <typename T1, typename T2>
        void operator() ( int i, T1 const& s, CORK::vector<T2> x, CORK::vector<T2> y ) const {
          throw exception::not_implemented( "there is no nonlinear matrix given for the requested value_type" ) ;
        }
      } ;

      struct solve_type {
        void operator() ( CORK::vector<Arg> x, bool is_new_shift ) const {
          throw exception::not_implemented( "there is no nonlinear matrix given for the requested value_type" ) ;
        }
      } ;

      struct funs_type {
        void operator() ( Arg const& val, CORK::vector<Arg> x ) const {
          throw exception::not_implemented( "there is no nonlinear matrix given for the requested value_type" ) ;
        }
      } ;

    public:
      empty()
      : basis_( funs_, 1 )
      , coefs_( 0, 0, multiply_add_, solve_ )
      {}

    public:
      auto matrix_polynomial() const {
       return CORK::matrix_valued_function::make_matrix_polynomial( basis_, coefs_ );
      }

    private:
      multiply_add_type                                                                 multiply_add_ ;
      solve_type                                                                        solve_ ;
      funs_type                                                                         funs_ ;
      basis::set_of_functions<Arg, funs_type>                                           basis_ ;
      coefficient_matrices::matrices_by_functions< Arg, multiply_add_type, solve_type > coefs_ ;
  } ; // empty()

} } // namespace CORK::matrix_valued_function

#endif
