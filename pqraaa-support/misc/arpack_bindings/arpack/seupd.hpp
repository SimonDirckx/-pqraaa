//
// Copyright (c) 2015
// Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
//

#ifndef ARPACK_SEUPD_HPP
#define ARPACK_SEUPD_HPP

#include <arpack/data.hpp>
#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <cassert>

//
// ARPACK computational routines
//

// Value-type variants of saupd
#define ARPACK_DSEUPD FORTRAN_ID( dseupd )
#define ARPACK_SSEUPD FORTRAN_ID( sseupd )
#define ARPACK_DSEUPD_C FORTRAN_ID( dseupd_c )
#define ARPACK_SSEUPD_C FORTRAN_ID( sseupd_c )

extern "C" {

void ARPACK_SSEUPD( fortran_bool_t* RVEC, char* HOWMNY, fortran_bool_t* SELECT, float* D, float* Z, int* LDZ
                  , float* SIGMA, char const* BMAT, int const* N, char const* WHICH, const int* NEV, float const* TOL
                  , float* RESID, int const* NCV, float* V, const int* LDV, int* IPARAM, int* IPNTR, float* WORKD
                  , float* WORKL, int const* LWORKL, int* INFO
                  ) ;

void ARPACK_SSEUPD_C( int* RVEC, char* HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, float* D, float* Z, int* LDZ
                  , float* SIGMA, char const* BMAT, int const* N, char const* WHICH, const int* NEV, float const* TOL
                  , float* RESID, int const* NCV, float* V, const int* LDV, int* IPARAM, int* IPNTR, float* WORKD
                  , float* WORKL, int const* LWORKL, int* INFO
                  ) ;

void ARPACK_DSEUPD( fortran_bool_t* RVEC, char* HOWMNY, fortran_bool_t* SELECT, double* D, double* Z, int* LDZ
                  , double* SIGMA, char const* BMAT, int const* N, char const* WHICH, const int* NEV, double const* TOL
                  , double* RESID, int const* NCV, double* V, const int* LDV, int* IPARAM, int* IPNTR, double* WORKD
                  , double* WORKL, int const* LWORKL, int* INFO
                  ) ;

void ARPACK_DSEUPD_C( int* RVEC, char* HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, double* D, double* Z, int* LDZ
                  , double* SIGMA, char const* BMAT, int const* N, char const* WHICH, const int* NEV, double const* TOL
                  , double* RESID, int const* NCV, double* V, const int* LDV, int* IPARAM, int* IPNTR, double* WORKD
                  , double* WORKL, int const* LWORKL, int* INFO
                  ) ;

}

// C++ bindings

namespace ARPACK {

  namespace detail {

    template <typename T>
    struct seupd_driver {
    } ;

    template <>
    struct seupd_driver< float > {
      static void apply( int RVEC, char HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, float* D, float* Z, int LDZ, float SIGMA
                       , char BMAT, int N, char const* WHICH, int NEV, float TOL
                       , float* RESID, int NCV, float* V, int LDV, int* IPARAM, int* IPNTR, float* WORKD
                       , float* WORKL, int LWORKL, int& INFO )
      {
        ARPACK_SSEUPD_C( &RVEC, &HOWMNY, SELECT_C, SELECT, D, Z, &LDZ, &SIGMA, &BMAT, &N, WHICH
                     , &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO
                     ) ;
      }
    } ;

    template <>
    struct seupd_driver< double > {
      static void apply( int RVEC, char HOWMNY, int* SELECT_C, fortran_bool_t* SELECT, double* D, double* Z, int LDZ, double SIGMA
                       , char BMAT, int N, char const* WHICH, int NEV, double TOL
                       , double* RESID, int NCV, double* V, int LDV, int* IPARAM, int* IPNTR, double* WORKD
                       , double* WORKL, int LWORKL, int& INFO )
      {
        ARPACK_DSEUPD_C( &RVEC, &HOWMNY, SELECT_C, SELECT, D, Z, &LDZ, &SIGMA, &BMAT, &N, WHICH
                     , &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, &INFO
                     ) ;
      }
    } ;

  } // detail

  // Binding for eigenvalues of a single matrix, no shift-and-invert.
  template <typename Select, typename SelectC, typename D, typename Z, typename Problem, typename Data>
  int seupd( bool RVEC, SelectC& select_c, Select& select, D& d, Z& z, Problem const& problem, Data& data ) {
    typedef typename Data::value_type value_type ;

    assert( boost::numeric::bindings::size(select)==data.ncv_ ) ;
    assert( boost::numeric::bindings::size(select_c)==data.ncv_ || boost::numeric::bindings::size(select_c)==0 ) ;

    assert( !RVEC || boost::numeric::bindings::size_column(z)==data.nev_ ) ;
    assert( boost::numeric::bindings::size(d)==data.nev_ ) ;

    char how_many ;
    if (boost::numeric::bindings::size(select_c)==0) {
      how_many = 'A' ;
    } else {
      how_many = 'S' ;
    }

    detail::seupd_driver<value_type>::apply( int(RVEC), how_many
                                           , boost::numeric::bindings::begin_value(select_c)
                                           , boost::numeric::bindings::begin_value(select)
                                           , boost::numeric::bindings::begin_value(d)
                                           , boost::numeric::bindings::begin_value(z)
                                           , boost::numeric::bindings::stride_major(z), problem.sigma()
                                           , data.bmat, data.n_, data.which_, data.nev_, data.tol_, data.resid_, data.ncv_, data.v_, data.ldv_
                                           , data.iparam_, data.ipntr_, data.workd_, data.workl_, data.lworkl_, data.info_
                                           ) ;
    return data.info_ ;
  } // seupd()

} // namespace ARPACK

#endif

