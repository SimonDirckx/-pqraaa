//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linear_solver_sparse_hpp
#define cork_linear_solver_sparse_hpp

#include <cork/coefficient_matrices/sparse.hpp>
//#include <cork/coefficient_matrices/any_glas.hpp>
#include <cork/linear_solver/linear_solver.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/sparse.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/sparse.hpp>
#include <boost/numeric/bindings/mumps/mumps_driver.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace linear_solver {

  template <typename ValueType >
  class sparse
  {
    public:
      typedef ValueType                        value_type ;
      typedef int                              size_type ;

    private:
      typedef glas2::coo<value_type, MUMPS_INT, std::ptrdiff_t, 1 >     coo_matrix_type ;
      typedef boost::numeric::bindings::mumps::mumps< coo_matrix_type > mumps_type ;

    public:
      sparse( int size )
      : factor_( std::make_shared<coo_matrix_type>(size, size) )
      , mumps_( std::make_shared<mumps_type>() )
      {
        //mumps_->sym = 2 ;
        mumps_->icntl[2] = mumps_->icntl[3] = 0 ; // Disable debug printout from MUMPS
        //mumps_->icntl[13] += 100000 ; // Extra fill-in, usually not needed
//        mumps_->icntl[10] = 1 ; // Condition number estimation
      }

    public:
      typedef sparse< ValueType > clone_type ;
      clone_type clone() const {
        return clone_type( factor_->num_rows() ) ;
      }

    public:
      template <typename Accumulate>
      void prepare_solve( Accumulate const& accumulate ) {
        factor_->reset( factor_->num_rows(), factor_->num_columns() ) ;
        accumulate( *factor_ ) ;

        ::boost::numeric::bindings::mumps::matrix_integer_data( *mumps_, *factor_ ) ;
        ::boost::numeric::bindings::mumps::matrix_value_data( *mumps_, *factor_ ) ;
        mumps_->job = 4 ;
        int info = ::boost::numeric::bindings::mumps::driver( *mumps_ ) ;
        if (info==-6 || info==-10) throw exception::linear_solver_failure() ;
        assert( info>=0 ) ;
      }

      template <typename W>
      void solve( W& w ) const {
        assert( w.size()==factor_->num_rows() ) ;
        auto w_r = reshape( w, w.size(), 1, glas2::column_major() ) ;
        ::boost::numeric::bindings::mumps::rhs_sol_value_data( *mumps_, w_r ) ;
        mumps_->job = 3 ;
#ifndef NDEBUG
        int info =
#endif
          ::boost::numeric::bindings::mumps::driver( *mumps_ ) ;
        assert( info==0 ) ;
//        std::cout << "Estimations of condition number of sparse linear system " << mumps_->rinfog[9] << " and " << mumps_->rinfog[10] << std::endl ;
      } // solve()

    private:
      mutable std::shared_ptr< coo_matrix_type > factor_ ;
      mutable std::shared_ptr< mumps_type >      mumps_ ;
  } ; // linear_solver


  template <typename T, typename Matrices>
  struct linear_solver_traits< T, coefficient_matrices::sparse<Matrices> > {
    typedef sparse< typename std::common_type<T, typename coefficient_matrices::sparse<Matrices>::value_type>::type > type ;
    static type apply( coefficient_matrices::sparse<Matrices> const& d ) {
      assert( d.num_rows()==d.num_columns() ) ;
      return type( d.num_rows() ) ;
    }
  } ;

/*
  template <typename T, typename Matrices>
  struct linear_solver_traits< T, coefficient_matrices::any_glas<Matrices>, typename std::enable_if< glas2::is<glas2::SparseMatrix, typename coefficient_matrices::any_glas<Matrices>::matrix_type>::value >::type > {
    typedef sparse< typename std::common_type<T, typename coefficient_matrices::sparse<Matrices>::value_type>::type > type ;

    template <typename J>
    static type apply( J const& d ) {
      assert( d.num_rows()==d.num_columns() ) ;
      return type( CORK::deref(d).num_rows() ) ;
    }
  } ;
*/


} } // namespace CORK::coefficient_matrices_linear_solver

#endif
