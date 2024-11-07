//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bindings_mumps_solver_hpp
#define glas2_bindings_mumps_solver_hpp

#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/coordinate.hpp>
#include <boost/numeric/bindings/mumps/mumps_driver.hpp>
#include <type_traits>

namespace glas2 {

  //
  // General class for linear solver
  // Must provide:
  //   linear_solver( M& )
  //   reset( M& )
  //   template <typename X> solve( X& x )
  //
  template <typename M>
  class mumps_solver
  {
    private:
      typedef boost::numeric::bindings::mumps::mumps< M > mumps_type ;
      static_assert( is<CoordinateSparseMatrix,M>::value, "mumps_solver: matrix M must be sparse" )
      static_assert( M::index_base==1, "mumps_solver: M::index_base != 1" )

    public:
      mumps_solver( M& m )
      {
        assert( m.num_rows()==m.num_columns() ) ;

        mumps_.icntl[2] = mumps_.icntl[3] = 0 ;
        mumps_.job = 1 ; ::boost::numeric::bindings::mumps::matrix_integer_data( mumps_, m ) ;
        ::boost::numeric::bindings::mumps::driver( mumps_ ) ;

        reset( m ) ;
      }

    public:
      int reset() {
        return ::boost::numeric::bindings::mumps::driver( mumps_ ) ;
      }

      int reset_structure( M& m ) {
        mumps_.job = 1 ; ::boost::numeric::bindings::mumps::matrix_integer_data( mumps_, m ) ;
        return reset() ;
      }

      int reset( M& m ) {
        mumps_.job = 2 ; ::boost::numeric::bindings::mumps::matrix_value_data( mumps_, m ) ;
        return reset() ;
      }

      template <typename X>
      int solve( X x ) const {
        mumps_.job = 3 ; rhs_sol_value_data( mumps_, x ) ;
        return ::boost::numeric::bindings::mumps::driver( mumps_ ) ;
      }

    private:
      mutable mumps_type mumps_ ;
  } ;

} // namespace glas2

#endif
