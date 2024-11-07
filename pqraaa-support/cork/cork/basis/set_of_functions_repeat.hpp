//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_set_of_functions_repeat_hpp
#define cork_basis_set_of_functions_repeat_hpp

#include <cork/basis/set_of_functions.hpp>
#include <cork/utility/value_type.hpp>
#include <cassert>
#include <functional>
#include <type_traits>

namespace CORK { namespace basis {

  template <typename ValueType, int Repeat, typename Function>
  class set_of_functions_repeat
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
      explicit set_of_functions_repeat( Function const& function, int num_terms )
      : function_( function )
      , num_terms_( num_terms )
      {} // set_of_functions

    public:
      template <typename FunctionValues>
      void evaluate( argument_value_type const& arg, FunctionValues values ) const {
        //static_assert( std::is_convertible< ValueType, argument_value_type >::value, "" ) ;
        static_assert( std::is_convertible< result_value_type, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms_*Repeat ) ;
        function_( arg, values(glas2::range(0,num_terms_)) ) ;
        for (int i=1; i<Repeat; ++i) values(glas2::range(i*num_terms_,i*num_terms_+num_terms_)) = arg * values(glas2::range(i*num_terms_-num_terms_,i*num_terms_)) ;
      } // evaluate()

    public:
      int num_terms() const { return num_terms_*Repeat ; }
      Function const& function() const { return function_ ; }

    public:
      auto naked_functions() const {
        return set_of_functions<ValueType,Function>( function_, num_terms_ ) ;
      }

    private:
      Function        function_ ;
      int             num_terms_ ;
  } ; // class set_of_functions_repeat

} } // namespace CORK::basis

#endif
