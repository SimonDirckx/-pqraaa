//  (C) Copyright Karl Meerbergen 2022.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linear_solver_user_defined_hpp
#define cork_linear_solver_user_defined_hpp

#include <cork/linear_solver/linear_solver.hpp>
#include <cork/coefficient_matrices/matrices_by_functions.hpp>
#include <cork/utility/ref.hpp>
#include <type_traits>
#include <cassert>
#include <typeinfo>

namespace CORK { namespace linear_solver {

  template <typename T, typename CoefficientMatrices>
  class user_defined
  {
    public:
      typedef typename CORK::deref_type<CoefficientMatrices>::type  coefficient_matrices_type ;
      typedef typename coefficient_matrices_type::size_type         size_type ;

    public:
      typedef T value_type ;

    public:
      explicit user_defined( CoefficientMatrices const& coefficient_matrices )
      : coefficient_matrices_( coefficient_matrices )
      , is_new_shift_( true )
      {}

    public:
      size_type size() const { return CORK::deref(coefficient_matrices_).num_rows() ; }

      // Old codes
    public:
      void prepare_solve( value_type const& shift ) {
        shift_ = shift ;
        is_new_shift_ = true ;
      } // prepare_solve()

      template <typename V>
      void solve( V v ) const {
        static_assert( std::is_same< typename std::decay<V>::type::value_type, value_type >::value, "V must have the same value_type as solver" ) ;
        CORK::deref(coefficient_matrices_).solve( shift_, v, is_new_shift_ ) ;
        is_new_shift_ = false ;
      } // solve()

    public:
      // 2020 codes
      template <typename V, typename Options>
      void solve( value_type& shift, V v, bool is_new_shift, Options const& ) const {
        static_assert( std::is_same< typename std::decay<V>::type::value_type, value_type >::value, "V must have the same value_type as solver" ) ;
        CORK::deref(coefficient_matrices_).solve( shift, v, is_new_shift_ ) ;
      } // solve()

    public:
      typedef user_defined clone_type ;
      clone_type clone() const {
        return *this ;
      }

    private:
      CoefficientMatrices  coefficient_matrices_ ;
      mutable bool         is_new_shift_ ;
      value_type           shift_ ;
  } ; // user_defined


  template <typename T, typename MultiplyAdd, typename Solve>
  struct linear_solver_traits< T, coefficient_matrices::matrices_by_functions< T, MultiplyAdd, Solve> > {
    typedef user_defined< T, coefficient_matrices::matrices_by_functions< T, MultiplyAdd, Solve> > type ;

    static type apply( coefficient_matrices::matrices_by_functions< T, MultiplyAdd, Solve> const& c ) {
      assert( c.num_rows()==c.num_columns() ) ;
      return type( c ) ;
    }
  } ;

} } // namespace CORK::linear_solver


#endif
