//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_user_defined_matrix_hpp
#define cork_coefficient_matrices_user_defined_matrix_hpp

#include <cork/utility/is_coefficient_matrix.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <cassert>
#include <type_traits>
#include <glas2/vector.hpp>

namespace CORK { namespace coefficient_matrices {

  template <typename Matrix>
  class user_defined_matrix
  {
    public:
      typedef typename std::decay<Matrix>::type matrix_type ;
      typedef typename matrix_type::value_type  value_type ;
      typedef typename matrix_type::size_type   size_type ;

    public:
      user_defined_matrix( Matrix m )
      : m_( m )
      {}

    public:
      size_type num_rows() const { return m_.size() ; }
      size_type num_columns() const { return m_.size() ; }

    public:
      Matrix const& coefficient_matrix() const { return m_ ; }

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
      void multiply_add( X const& x, W w ) const {
        static_assert( glas2::is<glas2::ContiguousDenseVector,W>::value, "W should be a contiguous vector" ) ;
        assert( x.size()==num_rows() ) ;
        assert( w.size()==num_rows() ) ;
        m_.multiply_add( multiply_add_traits<X>::scalar(x), x.size(), multiply_add_traits<X>::pointer(x), boost::numeric::bindings::begin_value(w) ) ;
      }

    private:
      Matrix m_ ;
  } ; // user_defined

  template <typename Matrix>
  user_defined_matrix<Matrix> make_user_defined_matrix( Matrix const& m ) { return user_defined_matrix<Matrix>( m ) ; }

  template <typename Matrix, typename EnableIf=void>
  struct is_user_defined_matrix
  : std::false_type
  {} ;

  template <typename Matrix>
  struct is_user_defined_matrix< user_defined_matrix<Matrix> >
  : std::true_type
  {} ;

} } // namespace CORK::coefficient_matrices


namespace CORK {

  template <typename Matrix>
  struct is_coefficient_matrix< coefficient_matrices::user_defined_matrix<Matrix> >
  : std::true_type
  {} ;

} // namespace CORK

#endif
