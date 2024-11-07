//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompany_glasing file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_norm_est_hpp
#define cork_coefficient_matrices_norm_est_hpp

#include <boost/numeric/bindings/lapack/computational/geqrf.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#include <cmath>

namespace CORK { namespace coefficient_matrices {

  template <typename T>
  class norm_est
  {
    public:
      //typedef typename std::decay<Matrices>::type matrices_type ;
      typedef T  value_type ;
      typedef decltype( std::abs(value_type()) )  real_type ;

    public:
      template <typename Matrices>
      norm_est( Matrices const& C )
      : R_( std::min<typename std::common_type<typename Matrices::grade_type, typename Matrices::size_type>::type>(C.num_rows(), C.num_matrices()), C.num_matrices() )
      //, n_( C.num_rows() )
      {
        assert( C.num_rows()==C.num_columns() ) ;
        glas2::vector< value_type > x( C.num_rows() ) ;
        glas2::matrix< value_type > Q( C.num_rows(), C.num_matrices() ) ;

        glas2::randomize( x ) ;
        glas2::fill( Q, 0.0 ) ;

        for (int i=0; i<C.num_matrices(); ++i) {

          C.multiply_add( i, x, Q(glas2::all(), i) ) ;
        }

        if (C.num_rows()<=R_.num_rows()) {
          R_ = Q ;
        } else {
          glas2::vector< value_type > tau( R_.num_rows() ) ;
#ifndef NDEBUG
          int info = 
#endif
            boost::numeric::bindings::lapack::geqrf( Q, tau ) ;
          assert( 0==info) ;
          fill( R_, 0.0 ) ;
          for (int i=0; i<C.num_matrices(); ++i) {
            R_( glas2::range(0,i+1), i ) = Q( glas2::range(0,i+1), i ) ;
          //R_ = Q( glas2::range(0,R_.num_rows()), glas2::all() ) ;
            //std::cout << "Matrix norm " << i << " = " << norm_2(R_( glas2::range(0,i+1), i ) ) << std::endl ;
          }
        }
      }

    protected:
      norm_est( int nrows, int ncols )
      : R_( nrows, ncols )
      {}

    public:
      real_type norm_coefficient( int i ) const {
        return norm_2( R_(glas2::all(), i ) ) ;
      } // norm_coefficient()

      template <typename Coefs>
      real_type norm( Coefs const& coefs ) const {
        glas2::vector< decltype(coefs(0)*R_(0,0)) > x( R_.num_rows() ) ;
        x = multiply( R_, coefs ) ;
        //return std::sqrt( real_type(n_) ) * norm_2( x ) ;
        return norm_2( x ) ;
      } // norm()

      template <typename Coefs>
      real_type upper_bound( Coefs const& coefs ) const {
        glas2::vector< decltype(std::abs(coefs(0))*std::abs(R_(0,0))) > x( R_.num_rows() ) ;
        x = multiply( abs(R_), abs(coefs) ) ;
        //return std::sqrt( real_type(n_) ) * norm_2( x ) ;
        return norm_2( x ) ;
      } // upper_bound()

    public:
      auto const& R() const { return R_ ; }

    protected:
      glas2::matrix<value_type> R_ ;
      //int                       n_ ;
  } ; // norm_est

} } // namespace CORK::coefficient_matrices

#endif
