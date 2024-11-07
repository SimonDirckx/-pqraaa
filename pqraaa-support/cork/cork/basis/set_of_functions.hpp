//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_set_of_functions_hpp
#define cork_basis_set_of_functions_hpp

#include <cork/concept/convertible.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/basis/iterator.hpp>
#include <cassert>
#include <functional>
#include <type_traits>

namespace CORK { namespace basis {

  template <typename ValueType, typename Function>
  class set_of_functions
  {
    public:
      typedef ValueType result_value_type ;
      typedef ValueType argument_value_type ;

      typedef int size_type ;

      template <typename T>
      using value_type_for = result_value_type ;

      template <typename T>
      using has_value_type_for = std::is_convertible<T, argument_value_type> ;

      template <typename Arg, typename T>
      using has_compatible_result_type = std::is_convertible<result_value_type, T> ;

      template <typename T>
      using has_compatible_argument_type = std::is_convertible<T, argument_value_type> ;

    public:
      explicit set_of_functions( Function const& function, int num_terms )
      : function_( function )
      , num_terms_( num_terms )
      {} // set_of_functions

    public:
      template <typename FunctionValues>
      void evaluate( argument_value_type const& arg, FunctionValues values ) const {
        //static_assert( std::is_convertible< ValueType, argument_value_type >::value, "" ) ;
        static_assert( std::is_convertible< result_value_type, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms_ ) ;
        function_( arg, values ) ;
      } // evaluate

    public:
      size_type num_terms() const { return num_terms_ ; }
      Function const& function() const { return function_ ; }

    private:
      Function        function_ ;
      int             num_terms_ ;
  } ; // class set_of_functions

} } // namespace CORK::basis

#endif
