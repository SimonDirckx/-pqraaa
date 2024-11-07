//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_lapack_order_schur_hpp
#define cork_lapack_order_schur_hpp

#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/data_order.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>

#include <cassert>
#include <string>

#define GLAS_LAPACK_STGSEN FORTRAN_ID( stgsen )
#define GLAS_LAPACK_DTGSEN FORTRAN_ID( dtgsen )
#define GLAS_LAPACK_CTGSEN FORTRAN_ID( ctgsen )
#define GLAS_LAPACK_ZTGSEN FORTRAN_ID( ztgsen )

#define GLAS_LAPACK_STRSEN FORTRAN_ID( strsen )
#define GLAS_LAPACK_DTRSEN FORTRAN_ID( dtrsen )
#define GLAS_LAPACK_CTRSEN FORTRAN_ID( ctrsen )
#define GLAS_LAPACK_ZTRSEN FORTRAN_ID( ztrsen )

extern "C" {
/*  void GLAS_LAPACK_DTGSEN ( int const* IJOB, fortran_bool_t const* WANTQ, fortran_bool_t const* WANTZ, fortran_bool_t const* SELECT, int const* N, double* A, int const* LDA
                          , double* B, int const* LDB, double* ALPHAR, double* ALPHAI, double* BETA, double* Q, int const* LDQ, double* Z, int const* LDZ
                          , int* M, double* PL, double* PR, double* DIF, double* WORK, int const* LWORK, int* IWORK, int const* LIWORK, int* INFO ) ;
  void GLAS_LAPACK_STGSEN ( int const* IJOB, fortran_bool_t const* WANTQ, fortran_bool_t const* WANTZ, fortran_bool_t const* SELECT, int const* N, float* A, int const* LDA
                          , float* B, int const* LDB, float* ALPHAR, float* ALPHAI, float* BETA, float* Q, int const* LDQ, float* Z, int const* LDZ
                          , int* M, float* PL, float* PR, float* DIF, float* WORK, int const* LWORK, int* IWORK, int const* LIWORK, int* INFO ) ;
*/
  void GLAS_LAPACK_DTRSEN( char const* JOB, char const* COMPQ, fortran_bool_t const* SELECT, int const* N, double* T, int const* LDT, double* Q, int const* LDQ
                         , double* WR, double* WI, int* M, double* S, double* SEP, double* WORK, int const* LWORK, int* IWORK, int const* LIWORK, int* INFO ) ;
  void GLAS_LAPACK_STRSEN( char const* JOB, char const* COMPQ, fortran_bool_t const* SELECT, int const* N, float* T, int const* LDT, float* Q, int const* LDQ
                         , float* WR, float* WI, int* M, float* S, float* SEP, double* WORK, int const* LWORK, int* IWORK, int const* LIWORK, int* INFO ) ;
//#ifdef GLAS_COMPLEX
  /*
  void GLAS_LAPACK_ZTGSEN ( int const* IJOB, fortran_bool_t const* WANTQ, fortran_bool_t const* WANTZ, fortran_bool_t const* SELECT, int const* N, std::complex<double>* A, int const* LDA
                          , std::complex<double>* B, int const* LDB, std::complex<double>* ALPHA, std::complex<double>* BETA, std::complex<double>* Q, int const* LDQ, std::complex<double>* Z, int const* LDZ
                          , int* M, double* PL, double* PR, double* DIF, std::complex<double>* WORK, int const* LWORK, int* IWORK, int const* LIWORK, int* INFO ) ;

  void GLAS_LAPACK_CTGSEN ( int const* IJOB, fortran_bool_t const* WANTQ, fortran_bool_t const* WANTZ, fortran_bool_t const* SELECT, int const* N, std::complex<float>* A, int const* LDA
                          , std::complex<float>* B, int const* LDB, std::complex<float>* ALPHA, std::complex<float>* BETA, std::complex<float>* Q, int const* LDQ, std::complex<float>* Z, int const* LDZ
                          , int* M, float* PL, float* PR, float* DIF, std::complex<float>* WORK, int const* LWORK, int* IWORK, int const* LIWORK, int* INFO ) ;

  void GLAS_LAPACK_ZTRSEN( char const* JOB, char const* COMPQ, fortran_bool_t const* SELECT, int const* N, std::complex<double>* T, int const* LDT, std::complex<double>* Q, int const* LDQ
                         , std::complex<double>* W, int* M, double* S, double* SEP, std::complex<double>* WORK, int const* LWORK, int* INFO ) ;
   */
//#endif

}

namespace CORK { namespace lapack { namespace detail {

  template <typename T>
  struct order_schur {} ;

  template <>
  struct order_schur<double> {
    template <typename A, typename B, typename X, typename Y, typename E, typename S>
    int operator() ( A& a, B& b, X& x, Y& y, E& e, S const& select, int& m ) const {
      int const ijob = 0 ;
      fortran_bool_t true_( true ) ;

      namespace bindings = boost::numeric::bindings ;
      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< B >::value, "order_schur: B must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;
      static_assert( bindings::is_column_major< Y >::value, "order_schur: Y must be column_major" ) ;

      int n = a.num_rows() ;
      assert( n==b.num_rows() ) ;
      assert( n==a.num_columns() ) ;
      assert( n==b.num_columns() ) ;
      assert( n==select.size() ) ;

      int lda = bindings::stride_major(a) ;
      int ldb = bindings::stride_major(b) ;
      int ldx = bindings::stride_major(x) ;
      int ldy = bindings::stride_major(y) ;
      double pl, pr, dif ;

      int const lwork = 4 * n + 16 ;
      glas2::vector< double > work( lwork ) ;

      glas2::vector< double > alphar( n ) ;
      glas2::vector< double > alphai( n ) ;
      glas2::vector< double > beta( n ) ;

      int const liwork = 1;
      int iwork ;

      int info ;

      GLAS_LAPACK_DTGSEN( &ijob, &true_, &true_, bindings::begin_value(select), &n, bindings::begin_value(a), &lda
                        , bindings::begin_value(b), &ldb, bindings::begin_value(alphar), bindings::begin_value(alphai)
                        , bindings::begin_value(beta), bindings::begin_value(y), &ldy, bindings::begin_value(x), &ldx
                        , &m, &pl, &pr, &dif, bindings::begin_value(work), &lwork, &iwork, &liwork, &info ) ;
      glas2::real(e) = alphar / beta ;
      glas2::imag(e) = alphai / beta ;

      return info ;
    }

    template <typename A, typename X, typename E, typename S>
    int operator() ( A& a, X& x, E& e, S const& select, int& m ) const {
      char const job = 'N' ;
      char const compq = 'V' ;

      namespace bindings = boost::numeric::bindings ;

      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;

      int n = a.num_rows() ;
      assert( n==a.num_columns() ) ;
      assert( n==select.size() ) ;

      int lda = bindings::stride_major(a) ;
      int ldx = bindings::stride_major(x) ;
      double s, sep ;

      int const lwork = n ;
      glas2::vector< double > work( lwork ) ;

      glas2::vector< double > wr( n ) ;
      glas2::vector< double > wi( n ) ;

      int const liwork = 1;
      int iwork ;

      int info ;

      GLAS_LAPACK_DTRSEN( &job, &compq, bindings::begin_value(select), &n, bindings::begin_value(a), &lda
                        , bindings::begin_value(x), &ldx
                        , bindings::begin_value(wr), bindings::begin_value(wi)
                        , &m, &s, &sep, bindings::begin_value(work), &lwork, &iwork, &liwork, &info ) ;

      glas2::real(e) = wr ; glas2::imag(e) = wi ;

      return info ;
    }
  } ; // order_schur

  template <>
  struct order_schur<float> {
    template <typename A, typename B, typename X, typename Y, typename E, typename S>
    int operator() ( A& a, B& b, X& x, Y& y, E& e, S const& select, int& m ) const {
      int const ijob = 0 ;
      fortran_bool_t true_( true ) ;

      namespace bindings = boost::numeric::bindings ;
      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< B >::value, "order_schur: B must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;
      static_assert( bindings::is_column_major< Y >::value, "order_schur: Y must be column_major" ) ;

      int n = a.num_rows() ;
      assert( n==b.num_rows() ) ;
      assert( n==a.num_columns() ) ;
      assert( n==b.num_columns() ) ;
      assert( n==select.size() ) ;

      int lda = bindings::stride_major(a) ;
      int ldb = bindings::stride_major(b) ;
      int ldx = bindings::stride_major(x) ;
      int ldy = bindings::stride_major(y) ;
      double pl, pr, dif ;

      int const lwork = 4 * n + 16 ;
      glas2::vector< float > work( lwork ) ;

      glas2::vector< float > alphar( n ) ;
      glas2::vector< float > alphai( n ) ;
      glas2::vector< float > beta( n ) ;

      int const liwork = 1;
      int iwork ;

      int info ;

      GLAS_LAPACK_STGSEN( &ijob, &true_, &true_, bindings::begin_value(select), &n, bindings::begin_value(a), &lda
                        , bindings::begin_value(b), &ldb, bindings::begin_value(alphar), bindings::begin_value(alphai)
                        , bindings::begin_value(beta), bindings::begin_value(y), &ldy, bindings::begin_value(x), &ldx
                        , &m, &pl, &pr, &dif, bindings::begin_value(work), &lwork, &iwork, &liwork, &info ) ;
      glas2::real( e ) = alphar / beta ;
      glas2::imag( e ) = alphai / beta ;

      return info ;
    }

    template <typename A, typename X, typename E, typename S>
    int operator() ( A& a, X& x, E& e, S const& select, int& m ) const {
      char const job = 'N' ;
      char const compq = 'V' ;

      namespace bindings = boost::numeric::bindings ;

      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;

      int n = a.num_rows(a) ;
      assert( n==a.num_columns(a) ) ;
      assert( n==select.size() ) ;

      int lda = bindings::stride_major(a) ;
      int ldx = bindings::stride_major(x) ;
      float s, sep ;

      int const lwork = n ;
      glas2::vector< float > work( lwork ) ;

      glas2::vector< float > wr( n ) ;
      glas2::vector< float > wi( n ) ;

      int const liwork = 1;
      int iwork ;

      int info ;

      GLAS_LAPACK_STRSEN( &job, &compq, bindings::begin_value(select), &n, bindings::begin_value(a), &lda
                        , bindings::begin_value(x), &ldx
                        , bindings::begin_value(wr), bindings::begin_value(wi)
                        , &m, &s, &sep, bindings::begin_value(work), &lwork, &iwork, &liwork, &info ) ;

      glas2::real(e) = wr ; glas2::imag(e) = wi ;

      return info ;
    }
  } ; // order_schur

//#ifdef GLAS_COMPLEX
  template <>
  struct order_schur< std::complex<double> > {
    template <typename A, typename B, typename X, typename Y, typename E, typename S>
    int operator() ( A& a, B& b, X& x, Y& y, E& e, S const& select, int& m ) const {
      int const ijob = 0 ;
      fortran_bool_t true_( true ) ;

      namespace bindings = boost::numeric::bindings ;
      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< B >::value, "order_schur: B must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;
      static_assert( bindings::is_column_major< Y >::value, "order_schur: Y must be column_major" ) ;

      int n = a.num_rows() ;
      assert( n==b.num_rows() ) ;
      assert( n==a.num_columns() ) ;
      assert( n==b.num_columns() ) ;
      assert( n==e.size() ) ;
      assert( n==select.size() ) ;
      assert( n==x.num_rows() ) ;
      assert( n==y.num_rows() ) ;
      assert( n==x.num_columns() ) ;
      assert( n==y.num_columns() ) ;

      int lda = bindings::stride_major(a) ;
      int ldb = bindings::stride_major(b) ;
      int ldx = bindings::stride_major(x) ;
      int ldy = bindings::stride_major(y) ;
      double pl, pr ;
      double dif[2];

      int const lwork = 1 ;
      std::complex<double> work ;

      glas2::vector< std::complex<double> > alpha( n ) ;
      glas2::vector< std::complex<double> > beta( n ) ;

      int const liwork = 1;
      int iwork ;

      int info ;

      GLAS_LAPACK_ZTGSEN( &ijob, &true_, &true_, bindings::begin_value(select), &n
                        , bindings::begin_value(a), &lda
                        , bindings::begin_value(b), &ldb, bindings::begin_value(alpha)
                        , bindings::begin_value(beta)
                        , bindings::begin_value(y), &ldy
                        , bindings::begin_value(x), &ldx
                        , &m, &pl, &pr, dif, &work, &lwork, &iwork, &liwork, &info ) ;
      e = alpha / beta ;

      return info ;
    }

    template <typename A, typename X, typename E, typename S>
    int operator() ( A& a, X& x, E& e, S const& select, int& m ) const {
      char const job = 'N' ;
      char const compq = 'V' ;

      namespace bindings = boost::numeric::bindings ;
      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;

      int n = a.num_rows() ;
      assert( n==a.num_columns() ) ;
      assert( n==select.size() ) ;

      int lda = bindings::stride_major(a) ;
      int ldx = bindings::stride_major(x) ;
      double s, sep ;

      int const lwork = 1 ;
      std::complex<double> work ;

      int info ;

      GLAS_LAPACK_ZTRSEN( &job, &compq, bindings::begin_value(select), &n, bindings::begin_value(a), &lda
                        , bindings::begin_value(x), &ldx
                        , bindings::begin_value(e)
                        , &m, &s, &sep, &work, &lwork, &info ) ;

      return info ;
    }
  } ; // order_schur

  template <>
  struct order_schur< std::complex<float> > {
    template <typename A, typename B, typename X, typename Y, typename E, typename S>
    int operator() ( A& a, B& b, X& x, Y& y, E& e, S const& select, int& m ) const {
      int const ijob = 0 ;
      fortran_bool_t true_( true ) ;

      namespace bindings = boost::numeric::bindings ;
      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< B >::value, "order_schur: B must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;
      static_assert( bindings::is_column_major< Y >::value, "order_schur: Y must be column_major" ) ;

      int n = a.num_rows() ;
      assert( n==b.num_rows() ) ;
      assert( n==a.num_columns() ) ;
      assert( n==b.num_columns() ) ;
      assert( n==select.size() ) ;

      int lda = bindings::stride_major(a) ;
      int ldb = bindings::stride_major(b) ;
      int ldx = bindings::stride_major(x) ;
      int ldy = bindings::stride_major(y) ;
      double pl, pr, dif ;

      int const lwork = 1 ;
      std::complex<float> work( lwork ) ;

      glas2::vector< std::complex<float> > alpha( n ) ;
      glas2::vector< std::complex<float> > beta( n ) ;

      int const liwork = 1;
      int iwork ;

      int info ;

      GLAS_LAPACK_CTGSEN( &ijob, &true_, &true_, bindings::begin_value(select), &n
                        , bindings::begin_value(a), &lda
                        , bindings::begin_value(b), &ldb, bindings::begin_value(alpha)
                        , bindings::begin_value(beta)
                        , bindings::begin_value(y), &ldy
                        , bindings::begin_value(x), &ldx
                        , &m, &pl, &pr, &dif, &work, &lwork, &iwork, &liwork, &info ) ;
      e = alpha / beta ;

      return info ;
    }

    template <typename A, typename X, typename E, typename S>
    int operator() ( A& a, X& x, E& e, S const& select, int& m ) const {
      char const job = 'N' ;
      char const compq = 'V' ;

      namespace bindings = boost::numeric::bindings ;
      static_assert( bindings::is_column_major< A >::value, "order_schur: A must be column_major" ) ;
      static_assert( bindings::is_column_major< X >::value, "order_schur: X must be column_major" ) ;

      int n = a.num_rows() ;
      assert( n==a.num_columns() ) ;
      assert( n==select.size() ) ;

      int lda = bindings::stride_major(a) ;
      int ldx = bindings::stride_major(x) ;
      double s, sep ;

      int const lwork = 1 ;
      std::complex<float> work ;

      int info ;

      GLAS_LAPACK_CTRSEN( &job, &compq, bindings::begin_value(select), &n, bindings::begin_value(a), &lda
                        , bindings::begin_value(x), &ldx
                        , bindings::begin_value(e)
                        , &m, &s, &sep, &work, &lwork, &info ) ;

      return info ;
    }
  } ; // order_schur
//#endif

} } } // namespace CORK::lapack::detail


namespace CORK { namespace lapack {

  template <typename A, typename B, typename X, typename Y, typename E, typename Keep>
  typename std::enable_if< glas2::is< glas2::DenseVector, Keep>::value, int >::type order_schur( A& a, B& b, X& x, Y& y, E& e, Keep const& keep, int& size ) {
    glas2::vector< fortran_bool_t > select( a.num_rows() ) ;
#ifndef NDEBUG
    for ( int i=0; i<keep.size(); ++i ) {
      assert( keep(i) < a.num_rows() ) ;
    }
#endif
    fill( select, false ) ;
    fill( select( keep ), true ) ;

    return detail::order_schur<typename A::value_type>() ( a, b, x, y, e, select, size ) ;
  } // order_schur()

  template <typename A, typename X, typename E, typename Keep>
  typename std::enable_if< glas2::is< glas2::DenseVector, Keep>::value, int >::type order_schur( A& a, X& x, E& e, Keep const& keep, int& size ) {
    glas2::vector< fortran_bool_t > select( a.num_rows() ) ;
#ifndef NDEBUG
    for ( int i=0; i<keep.size(); ++i ) {
      assert( keep(i) < a.num_rows() ) ;
    }
#endif
    fill( select, false ) ;
    fill( select( keep ), true ) ;

    return detail::order_schur<typename A::value_type>() ( a, x, e, select, size ) ;
  } // order_schur()

  template <typename A, typename B, typename X, typename Y, typename E>
  int order_schur( A& a, B& b, X& x, Y& y, E& e, std::string const& todo, int& size ) {
    glas2::vector< fortran_bool_t > select( a.num_rows() ) ;
    fill( select, false ) ;

    if (todo=="real") {
      for ( int i=0; i<e.size(); ++i ) {
        if ( glas2::imag(e(i))==0.0 ) select(i) = fortran_bool_t(true) ;
      }
    } else {
      std::cerr << std::string("String ") + todo + " not known to ORDER_SCHUR"  << std::endl ;
      return -1 ;
    }

    return detail::order_schur<typename A::value_type>() ( a, b, x, y, e, select, size ) ;
  }

  template <typename A, typename X, typename E>
  int order_schur( A& a, X& x, E& e, std::string const& todo, int& size ) {
    glas2::vector< fortran_bool_t > select( a.num_rows() ) ;
    fill( select, false ) ;

    if (todo=="real") {
      for ( int i=0; i<e.size(); ++i ) {
        if ( glas2::imag(e(i))==0.0 ) select(i) = fortran_bool_t(true) ;
      }
    } else {
      std::cerr << std::string("String ") + todo + " not known to ORDER_SCHUR"  << std::endl ;
      return -1 ;
    }

    return detail::order_schur<typename A::value_type>() ( a, x, e, select, size ) ;
  }

} } // namespace CORK::lapack

#endif
