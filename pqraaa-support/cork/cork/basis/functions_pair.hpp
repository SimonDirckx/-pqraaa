//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_functions_pair_hpp
#define cork_basis_functions_pair_hpp

#include <cork/exception/not_implemented.hpp>
#include <cork/basis/iterator.hpp>
#include <cassert>
#include <functional>
#include <type_traits>

namespace CORK { namespace basis {

  template <typename Functions1, typename Functions2>
  class functions_pair
  {
    public:
      typedef Functions1 functions_1_type ;
      typedef Functions2 functions_2_type ;
      static_assert( !std::is_same< typename functions_1_type::argument_value_type, typename functions_2_type::argument_value_type >::value
                   , "CORK:: scalar_functions_pair: argument_value_types must be different" ) ;

      typedef int size_type ;

      template <typename T>
      using value_type_for = typename std::conditional< typename functions_1_type::template has_value_type_for<T>::value
                                                      , typename functions_1_type::template value_type_for<T>::value
                                                      , typename functions_2_type::template value_type_for<T>::value
                                                      >::type ;

      template <typename T>
      using has_value_type_for = std::bool_constant< std::conditional< typename functions_1_type::template has_value_type_for<T>::value
                                                  && std::conditional< typename functions_2_type::template has_value_type_for<T>::value
                                                   > ;

      template <typename Arg, typename T>
      using has_compatible_result_type = std::is_convertible< typename value_type_for<Arg>::type, T> ;

      template <typename T>
      using has_compatible_argument_type = has_value_type_for<T> ;

    public:
      explicit functions_pair( functions_1_type const& functions_1, functions_2_type const& functions_2
      : functions_1_( functions_1 )
      , functions_2_( functions_2 )
      {} // functions_pair

    private:
      template <typename T>
      struct evaluate_traits {
        template <typename Arg, typename FunctionValues>
        static void apply( Arg const& arg, FunctionValues& values ) {
          throw exception::not_implemented("CORK is requesting scalar_function evaluation for another value type, which is not provided by you. Sorry") ;
        }
      } ;

      template <typename functions_1_type::argument_value_type>
      struct evaluate_traits {
        template <typename Arg, typename FunctionValues>
        static void apply( Arg const& arg, FunctionValues& values ) {
          functions_1_.evaluate( arg, values ) ;
        }
      } ;

      template <typename functions_2_type::argument_value_type>
      struct evaluate_traits {
        template <typename Arg, typename FunctionValues>
        static void apply( Arg const& arg, FunctionValues& values ) {
          functions_2_.evaluate( arg, values ) ;
        }
      } ;

    public:
      template <typename ValueType, typename FunctionValues>
      void evaluate( ValueType const& arg, FunctionValues values ) const {
        assert( values.size() == num_terms_ ) ;
        function_( arg, values ) ;
      } // evaluate

    public:
      int const& num_terms() const {
        assert( functions_1_.num_terms() == functions_2_.num_terms() ) ;
        return functions_1_num_terms_ ;
      }

    private:
      Function const& function_ ;
      int             num_terms_ ;
  } ; // class functions_pair

} } // namespace CORK::basis

#endif
