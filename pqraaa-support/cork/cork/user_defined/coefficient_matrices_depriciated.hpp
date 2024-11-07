//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_user_defined_coefficient_matrices_hpp
#define cork_user_defined_coefficient_matrices_hpp

#include <cork/coefficient_matrices/coefficient_matrices.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <cassert>
#include <type_traits>
#include <glas2/vector.hpp>

namespace CORK { namespace user_defined {

  template <typename ValueType, typename MultiplyAdd, typename Solve>
  class coefficient_matrices
  {
    public:
      typedef int       size_type ;
      typedef int       grade_type ;
      typedef ValueType value_type ;

    public:
      coefficient_matrices( size_type n, size_type d, MultiplyAdd const& multiply_add, Solve& solve )
      : n_( n )
      , d_( d )
      , multiply_add_( multiply_add )
      , solve_( solve )
      {}

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
        static auto pointer( V const& v ) { return boost::numeric::bindings::begin_value(v) ; }
      } ;

      template <typename S, typename V, typename Op>
      struct multiply_add_traits< glas2::binary_operation<S,V,Op>, typename std::enable_if< glas2::is<glas2::ContiguousDenseVector,V>::value >::type > {
        static_assert( glas2::is<glas2::ContiguousDenseVector,V>::value, "V should be a contiguous vector" ) ;
        static auto scalar( glas2::binary_operation<S,V,Op> const& v ) { return v.scalar() ; }
        static auto pointer( glas2::binary_operation<S,V,Op> const& v ) { return boost::numeric::bindings::begin_value(v.vector()) ; }
      } ;

      template <typename S, typename V, typename Op>
      struct multiply_add_traits< glas2::binary_operation<S,V,Op>, typename std::enable_if< !glas2::is<glas2::ContiguousDenseVector,V>::value >::type > {
        static auto scalar( glas2::binary_operation<S,V,Op> const& v ) { return v.scalar()*multiply_add_traits<V>::scalar(v.vector()) ; }
        static auto pointer( glas2::binary_operation<S,V,Op> const& v ) { return multiply_add_traits<V>::pointer(v.vector()) ; }
      } ;

    public:
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        static_assert( glas2::is<glas2::ContiguousDenseVector,W>::value, "W should be a contiguous vector" ) ;
        assert( x.size()==num_rows() ) ;
        assert( w.size()==num_rows() ) ;
        assert( i>=0 && i<num_matrices() ) ;
        multiply_add_( i, multiply_add_traits<X>::scalar(x), x.size(), multiply_add_traits<X>::pointer(x), boost::numeric::bindings::begin_value(w) ) ;
      } // multiply_add()

      template <typename Shift, typename X>
      void solve( Shift const& shift, X x, bool is_new_shift ) {
        static_assert( glas2::is<glas2::ContiguousDenseVector,X>::value, "X should be a contiguous vector" ) ;
        assert( x.size()==num_rows() ) ;
        solve_( shift, x.size(), boost::numeric::bindings::begin_value(x), is_new_shift ) ;
      } // solve()

    private:
      size_type          n_ ;
      size_type          d_ ;
      MultiplyAdd const& multiply_add_ ;
      Solve      &       solve_ ;
  } ; // coefficient_matrices


} } // namespace CORK::user_defined

#endif
