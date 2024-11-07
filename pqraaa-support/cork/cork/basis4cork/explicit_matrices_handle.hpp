//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_explicit_matrices_handle_hpp
#define cork_basis4cork_explicit_matrices_handle_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getrs.hpp>
//#include <boost/numeric/bindings/lapack/computational/gecon.hpp>
#include <cassert>
#include <string>
#include <stdexcept>
#include <functional>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename ValueType, typename T>
  class explicit_matrices_handle
  {
    public:
      typedef ValueType                  value_type ;

    private:
      typedef typename std::common_type< T, ValueType >::type  inner_value_type ;
      typedef glas2::shared_matrix< T >                        matrix_type ;
      typedef glas2::shared_matrix< inner_value_type >         matrix_temp_type ;
      typedef glas2::shared_vector< int >                      int_vector_type ;

    public:
      typedef typename matrix_type::size_type size_type ;
    
    public:
      explicit explicit_matrices_handle( matrix_type const& M, matrix_type const& N, std::function<value_type(value_type)> const& phi0 )
      : M_( M )
      , N_( N )
      , lu_( M.num_rows(), M.num_rows() )
      , pivots_( M.num_rows() )
      , phi_0_( phi0 )
      {}

      explicit explicit_matrices_handle( matrix_type const& M, matrix_type const& N )
      : M_( M )
      , N_( N )
      , lu_( M.num_rows(), M.num_rows() )
      , pivots_( M.num_rows() )
      , phi_0_( []( value_type const& ) { return 1.; } )
      {}

    protected:
      auto const& M() const { return M_ ; }
      auto const& N() const { return N_ ; }

      typename matrix_type::base_type M() { return M_ ; }
      typename matrix_type::base_type N() { return N_ ; }

    public:
      size_type size() const { return M_.num_columns() ; }

      value_type phi_0() const { return phi_0_( s_ ) ; }

      // Basis4CORK
      void shift( ValueType s ) {
        s_ = s ;

        fill(lu_,0.0) ;
        lu_ = M_(glas2::all(), glas2::range_from_end(1,0)) - s_ * N_(glas2::all(), glas2::range_from_end(1,0)) ;
        //auto anorm = norm_1( lu_ ) ;

        int info = boost::numeric::bindings::lapack::getrf( lu_, pivots_ ) ;
        assert( info>=0 ) ;
        if (info!=0) throw std::runtime_error( "CORK::basis4CORK::explicit_matrices: shift failed") ;
        /*decltype(anorm) rcond ;
        info = boost::numeric::bindings::lapack::gecon( '1', lu_, anorm, rcond ) ;
        std::cout << "M-s*N condition number estimation " << 1./rcond << std::endl ;
        */
      }

      // Basis4CORK
      ValueType shift() const { return s_ ; }

      // Basis4CORK
      template <typename Z0, typename ZZ, typename Backend=glas2::default_backend>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z, Backend const& backend=glas2::default_backend() ) const {
        assert( Z.num_rows()==M_.num_rows() ) ;
        fill( Z, 0.0 ) ;
        Z = outer_prod( M_(glas2::all(),0), z0 ) - s_ * outer_prod( N_(glas2::all(),0), z0 ) ;
      } // lower_solve_right_hand_side()

      // Basis4CORK
      template <typename ZM, typename Backend=glas2::default_backend>
      void solve( ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assert( Z.num_rows()==lu_.num_rows() ) ;
        glas2::matrix< typename ZM::value_type > Zc(Z.num_rows(), Z.num_columns() ) ;
        Zc = Z ;
        boost::numeric::bindings::lapack::getrs( lu_, pivots_, Zc ) ;
        Z = Zc ;
      } // solve()

    private:
      matrix_type const&                    M_ ;
      matrix_type const&                    N_ ;
      matrix_temp_type                      lu_ ;
      int_vector_type                       pivots_ ;
      value_type                            s_ ;
      std::function<value_type(value_type)> phi_0_ ;
  } ; // explicit_matrices_handle

} } // namespace CORK::basis4cork

#endif
