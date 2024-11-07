//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_tuple_of_functions_hpp
#define cork_approximation_tuple_of_functions_hpp

#include <cassert>
#include <type_traits>
#include <tuple>

namespace CORK { namespace approximation {

  namespace detail {

    template <int I, typename Tuple, typename Arg, typename FunctionValues>
    typename std::enable_if< (std::tuple_size<Tuple>::value == I) >::type evaluate_1( Tuple const& tuple, Arg const& arg, FunctionValues values ) {
    }

    template <int I, typename Tuple, typename Arg, typename FunctionValues>
    typename std::enable_if< (std::tuple_size<Tuple>::value > I) >::type evaluate_1( Tuple const& tuple, Arg const& arg, FunctionValues values ) {
      values(I) = std::get<I>(tuple)( arg ) ;
      evaluate_1<I+1>( tuple, arg, values ) ;
    }
  } // namespace detail


  template <typename Tuple, typename T=typename std::tuple_element<0,Tuple>::type::result_type>
  class tuple_of_functions
  {
    public:
      typedef typename std::tuple_size<Tuple>::value_type size_type ;
      typedef T value_type ;

    public:
      explicit tuple_of_functions( Tuple const& tuple )
      : tuple_( tuple )
      {}

    public:
      size_type num_terms() const {
        return std::tuple_size<Tuple>::value ;
      } // num_terms()

    public:
      template <typename ValueType, typename FunctionValues>
      void evaluate( ValueType const& arg, FunctionValues values ) const {
        assert( values.size() == num_terms() ) ;
        detail::evaluate_1<0>( tuple_, arg, values ) ;
      } // evaluate

    private:
      Tuple tuple_ ;

  } ; // class tuple_of_functions

  template <typename Tuple>
  tuple_of_functions<Tuple> make_tuple_of_functions( Tuple const& tuple ) {
    return tuple_of_functions<Tuple>( tuple ) ;
  }

} } // namespace CORK::approximation

#endif
