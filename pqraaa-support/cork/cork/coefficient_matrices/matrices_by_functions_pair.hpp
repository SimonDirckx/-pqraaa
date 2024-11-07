//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_matrices_by_functions_pair_hpp
#define cork_coefficient_matrices_matrices_by_functions_pair_hpp

#include <cork/exception/not_implemented.hpp>
#include <cork/coefficient_matrices/matrices_by_functions.hpp>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  template <typename Matrices1, typename Matrices2>
  class matrices_by_functions_pair
  {
    public:
      typedef int       size_type ;
      typedef int       grade_type ;
      static_assert( !std::is_same< typename Matrices1::value_type, typename Matrices2::value_type >::value ) ;

      template <typename T>
      using value_type_for = typename std::conditional< Matrices1::template has_value_type_for<T>::value
                                                      , Matrices1::template value_type_for<T>::type
                                                      , Matrices2::template value_type_for<T>::type
                                                      >::type ;

      template <typename T>
      using has_value_type_for = std::bool_constant< Matrices1::template has_value_type_for<T>::value || Matrices2::template ha_value_type_for<T>::value > ;

    public:
      matrices_by_functions( Matrices1 const& matrices_1, Matrices2 const& matrices2 )
      : matrices_1_( matrices1 )
      , matrices_2_( matrices2 )
      {}

    public:
      size_type num_rows() const {
        assert( matrices_1_.num_rows()==matrices_2_.num_rows() ) ;
        return matrices_1_.num_rows() ;
      }
      size_type num_columns() const {
        assert( matrices_1_.num_columns()==matrices_2_.num_columns() ) ;
        return matrices_1_.num_columns() ;
      }

      grade_type num_matrices() const {
        assert( matrices_1_.num_matrices()==matrices_2_.num_matrices() ) ;
        return matrices_1_.num_matrices() ;
      }

    private:
      template <typename T>
      struct multiply_add_traits {
        template <typename X, typename W>
        static void apply( size_type i, X const& x, W& w ) {
          throw exception::not_implemented("CORK is requesting 'multiply_add' for another value type, which is not provided by you. Sorry") ;
        }
      } ;

      template <>
      struct multiply_add_traits< typename Matrices1::value_type > {
        template <typename X, typename W>
        static void apply( size_type i, X const& x, W& w ) {
          matrices_1_.multiply_add( i, x, w ) ;
        }
      } ;

      template <>
      struct multiply_add_traits< typename Matrices2::value_type > {
        template <typename X, typename W>
        static void apply( size_type i, X const& x, W& w ) {
          matrices_2_.multiply_add( i, x, w ) ;
        }
      } ;

    public:
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        static_assert( has_value_type_for<typename X::value_type>::value, "CORK::multiply_add: value_type is not compatible with this computation" ) ;
        static_assert( has_value_type_for<typename W::value_type>::value, "CORK::multiply_add: value_type is not compatible with this computation" ) ;
        static_assert( glas2::is<glas2::ContiguousDenseVector,W>::value, "W should be a contiguous vector" ) ;
        multiply_add_traits<typename W::value_type>::apply( i x, w ) ;
      } // multiply_add()

    private:
      template <typename T>
      struct solve_traits {
        template <typename Shift, typename X>
        static void apply( Shift const& shift, X x, bool is_new_shift ) {
          throw exception::not_implemented("CORK is requesting 'solve' for another value type, which is not provided by you. Sorry") ;
        }
      } ;

      template <>
      struct solve_traits< typename Matrices1::value_type > {
        template <typename Shift, typename X>
        static void apply( Shift const& shift, X x, bool is_new_shift ) {
          matrices_1_.solve( shift, x, is_new_shift ) ;
        }
      } ;

      template <>
      struct solve_traits< typename Matrices2::value_type > {
        template <typename Shift, typename X>
        static void apply( Shift const& shift, X x, bool is_new_shift ) {
          matrices_2_.solve( shift, x, is_new_shift ) ;
        }
      } ;

    public:
      template <typename Shift, typename X>
      void solve( Shift const& shift, X x, bool is_new_shift ) {
        static_assert( has_value_type_for<typename X::value_type>::value, "CORK::multiply_add: value_type is not compatible with this computation" ) ;
        static_assert( glas2::is<glas2::ContiguousDenseVector,X>::value, "X should be a contiguous vector" ) ;
        assert( x.size()==num_rows() ) ;
        solve_traits<typename X::value_type>::apply( shift, x, is_new_shift ) ;
      } // solve()

    private:
      Matrices1 const& matrices_1_ ;
      Matrices2 const& matrices_2_ ;
  } ; // coefficient_matrices


} } // namespace CORK::coefficient_matrices

#endif
