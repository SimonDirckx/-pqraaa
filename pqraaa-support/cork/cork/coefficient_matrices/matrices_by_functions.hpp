//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_matrices_by_functions_hpp
#define cork_coefficient_matrices_matrices_by_functions_hpp

#include <cork/coefficient_matrices/coefficient_matrices.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <cork/utility/ref.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <cassert>
#include <iostream>
#include <type_traits>
#include <glas2/vector.hpp>

namespace CORK { namespace coefficient_matrices {

  template <typename ValueType, typename MultiplyAdd, typename Solve>
  class matrices_by_functions
  {
    public:
      typedef int       size_type ;
      typedef int       grade_type ;
      typedef ValueType value_type ;

      template <typename T>
      using value_type_for = value_type ;

      template <typename T>
      using has_value_type_for = std::is_same< value_type, T > ;

    public:
      matrices_by_functions( size_type n, size_type d, MultiplyAdd multiply_add, Solve solve )
      : n_( n )
      , d_( d )
      , multiply_add_( multiply_add )
      , solve_( solve )
      {
//#ifndef NDEBUG
//        std::cout << "CORK::nonlinear_matrix: arguments are copied" << std::endl ;
//#endif
      }

    public:
      size_type num_rows() const { return n_ ; }
      size_type num_columns() const { return n_ ; }

      grade_type num_matrices() const { return d_+1 ; }

    private:
      template <typename Expression, typename EnableIf=void>
      struct multiply_add_traits {
      } ;

      template <typename V>
      struct multiply_add_traits< V, typename std::enable_if< glas2::is<glas2::ContiguousDenseVector,V>::value >::type > {
        static auto scalar( V const& v ) { return 1.0 ; }
        static auto vector( V const& v ) { return v; }
      } ;

      template <typename S, typename V, typename Op>
      struct multiply_add_traits< glas2::binary_operation<S,V,Op>, typename std::enable_if< glas2::is<glas2::ContiguousDenseVector,V>::value >::type > {
        static_assert( glas2::is<glas2::ContiguousDenseVector,V>::value, "V should be a contiguous vector" ) ;
        static auto scalar( glas2::binary_operation<S,V,Op> const& v ) { return v.scalar() ; }
        static auto vector( glas2::binary_operation<S,V,Op> const& v ) { return v.vector() ;}
      } ;

/*      template <typename S, typename V, typename Op>
      struct multiply_add_traits< glas2::binary_operation<S,V,Op>, typename std::enable_if< !glas2::is<glas2::ContiguousDenseVector,V>::value >::type > {
        static auto scalar( glas2::binary_operation<S,V,Op> const& v ) { return v.scalar()*multiply_add_traits<V>::scalar(v.vector()) ; }
        static auto vector( glas2::binary_operation<S,V,Op> const& v ) { return v.vector() ;}
      } ;*/

    public:
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        static_assert( has_value_type_for<typename X::value_type>::value, "CORK::multiply_add: value_type is not compatible with this computation" ) ;
        static_assert( has_value_type_for<typename W::value_type>::value, "CORK::multiply_add: value_type is not compatible with this computation" ) ;
        static_assert( glas2::is<glas2::ContiguousDenseVector,W>::value, "W should be a contiguous vector" ) ;
        assert( x.size()==num_rows() ) ;
        assert( w.size()==num_rows() ) ;
        assert( i>=0 && i<num_matrices() ) ;
        auto x_ = multiply_add_traits<X>::vector(x) ;
        auto scalar_ = multiply_add_traits<X>::scalar(x) ;
        CORK::deref(multiply_add_)( i, scalar_, x_, w ) ;
      } // multiply_add()

      template <typename Shift, typename X>
      void solve( Shift const& shift, X x, bool is_new_shift ) const {
        static_assert( has_value_type_for<typename X::value_type>::value, "CORK::multiply_add: value_type is not compatible with this computation" ) ;
        static_assert( glas2::is<glas2::ContiguousDenseVector,X>::value, "X should be a contiguous vector" ) ;
        assert( x.size()==num_rows() ) ;
        CORK::deref(solve_)( shift, x, is_new_shift ) ;
      } // solve()

    public:
      template <typename Coefs, typename Matrix>
      typename std::enable_if< glas2::is<glas2::DenseMatrix, Matrix>::value>::type accumulate( Coefs const& coefs, Matrix& A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
        assert( coefs.size()==num_matrices() ) ;
        assert( A.num_rows() == n_ ) ;
        assert( A.num_columns() == n_ ) ;

        glas2::vector< typename std::common_type<typename Coefs::value_type, typename Matrix::value_type>::type > temp( n_ ) ;
        fill(temp, 0.0) ;

        for ( int j=0; j<n_; ++j ) {
          temp( j ) = 1.0 ;
          for ( int i=0; i<coefs.size(); ++i ) {
             CORK::deref(multiply_add_)( i, coefs(i), temp, A( glas2::all(), j ) ) ;
          }
          temp( j ) = 0.0 ;
        }
      } // accumulate()

    private:
      size_type    n_ ;
      size_type    d_ ;
      MultiplyAdd  multiply_add_ ;
      Solve        solve_ ;
  } ; // matrices_by_functions


  template <typename M>
  struct is_matrices_by_functions
  : std::false_type
  {} ;


  template <typename ValueType, typename MultiplyAdd, typename Solve>
  struct is_matrices_by_functions< matrices_by_functions<ValueType,MultiplyAdd,Solve> >
  : std::true_type
  {} ;


} } // namespace CORK::coefficient_matrices

#endif
