//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linear_solver_dense_hpp
#define cork_linear_solver_dense_hpp

#include <cork/coefficient_matrices/dense.hpp>
//#include <cork/coefficient_matrices/any_glas.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <cork/utility/ref.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getrs.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include <cassert>
#include <type_traits>
#include <typeinfo>

namespace CORK { namespace linear_solver {


  template <typename ValueType >
  class dense
  {
    public:
      typedef ValueType                        value_type ;
      typedef int                              size_type ;
      typedef glas2::shared_matrix<value_type> dense_matrix_type ;

    public:
      dense( int n )
      : factor_( n, n )
      , pivots_( n )
      {}

    public:
      template <typename Accumulate>
      void prepare_solve( Accumulate const& accumulate ) {
        fill( factor_, 0.0 ) ;
        accumulate( factor_ ) ;
        int info = boost::numeric::bindings::lapack::getrf( factor_, pivots_ ) ;
        assert( info>=0 ) ;
        if (info>0) throw exception::linear_solver_failure() ;
      } // prepare_solve()

      template <typename W>
      void solve( W& w ) const {
        static_assert( std::is_same<value_type, typename W::value_type>::value, "" ) ;
        assert( w.size()==factor_.num_rows() ) ;

#ifndef NDEBUG
        int info =
#endif
        boost::numeric::bindings::lapack::getrs( factor_, pivots_, w ) ;
        assert( info==0 ) ;
      } // solve()

    private:
      dense_matrix_type   factor_ ;
      std::vector<int>    pivots_ ;
  } ; // class dense


  template <typename T, typename Matrices>
  struct linear_solver_traits< T, coefficient_matrices::dense<Matrices> > {
    typedef dense< typename std::common_type<T, typename coefficient_matrices::dense<Matrices>::value_type>::type > type ;

    static type apply( coefficient_matrices::dense<Matrices> const& d ) {
      assert( d.num_rows()==d.num_columns() ) ;
      return type( d.num_rows() ) ;
    }
  } ;

  /* is defined in linear_solver::any_glas
  template <typename T, typename Matrices>
  struct linear_solver_traits< T, coefficient_matrices::any_glas<Matrices>, typename std::enable_if< glas2::is<glas2::DenseMatrix, typename coefficient_matrices::any_glas<Matrices>::matrix_type>::value >::type > {
    typedef dense< typename std::common_type<T, typename coefficient_matrices::any_glas<Matrices>::value_type>::type > type ;

    static type apply( coefficient_matrices::any_glas<Matrices> const& d ) {
      assert( d.num_rows()==d.num_columns() ) ;
      return type( d.num_rows() ) ;
    }
  } ;
  */

} } // namespace CORK::linear_solver

#endif
