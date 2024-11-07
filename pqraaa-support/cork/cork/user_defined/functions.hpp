//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_user_defined_functions_hpp
#define cork_user_defined_functions_hpp

#include <cork/utility/value_type_for.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/basis/iterator.hpp>
#include <cassert>
#include <type_traits>
#include <vector>
#include <functional>

namespace CORK { namespace user_defined {

  template <typename Sequence>
  class functions
  {
    public:
      typedef typename std::decay< Sequence >::type sequence_type ;
      typedef typename sequence_type::size_type     size_type ;

    public:
      explicit functions( Sequence sequence )
      : sequence_( sequence )
      {} // functions

    public:
      size_type num_terms() const {
        return sequence_.size() ;
      }

    public:
      template <typename ValueType, typename FunctionValues>
      void evaluate( ValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< ValueType, typename sequence_type::value_type::result_type >::value, "" ) ;
        static_assert( std::is_convertible< typename value_type<functions<Sequence>>::type, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;
        auto it = sequence_.begin() ;
        for (typename std::decay<FunctionValues>::type::size_type i=0; i<values.size(); ++i) {
          values(i) = (*it)( arg ) ;
          ++it ;
        }
      } // evaluate

    public:
      Sequence const& sequence() const { return sequence_ ; }

    private:
      Sequence sequence_ ;
  } ; // class functions

} } // namespace CORK::user_defined


namespace CORK {

  template <typename Sequence>
  struct value_type< user_defined::functions<Sequence> >
  : std::invoke_result< typename std::decay< Sequence >::type::value_type >
  //: std::result_of< typename std::decay< Sequence >::type::value_type >
  {
 //   typedef typename std::decay< Sequence >::type::value_type::result_type type ;
  } ;

  template <typename T, typename Sequence>
  struct value_type_for< T, user_defined::functions<Sequence> >
  : value_type< user_defined::functions<Sequence> >
  {} ;

} // namespace CORK

namespace std {
 
  template <typename T, typename S>
  struct invoke_result< std::function<T(S)> > {
    typedef T type ;
  } ;
}

#endif
