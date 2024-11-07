//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_functions_hpp
#define cork_basis_functions_hpp

#include <cork/utility/value_type_for.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/concept/convertible.hpp>
#include <cork/basis/iterator.hpp>
#include <cassert>
#include <functional>
#include <type_traits>

namespace CORK { namespace basis {

  template <typename Sequence>
  class functions
  {
    public:
      typedef typename std::decay< Sequence >::type sequence_type ;
      typedef typename sequence_type::size_type     size_type ;

      typedef typename std::invoke_result< typename std::decay< Sequence >::type::value_type >::type value_type ;

      template <typename T>
      using value_type_for = decltype( typename Sequence::value_type()( T() ) ) ;

      template <typename T>
      using has_value_type_for = CORK::is_convertible<T, value_type> ;

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
  //      static_assert( std::is_convertible< ValueType, typename sequence_type::value_type::result_type >::value, "" ) ;
 //       static_assert( std::is_convertible< typename value_type<functions<Sequence>>::type, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        static_assert( std::is_convertible< value_type, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( (long int)(values.size()) == (long int)(num_terms()) ) ;
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

  template <typename Sequence>
  decltype(auto) make_sequence_of_functions( Sequence const& sequence ) { return functions<Sequence>( sequence ) ; }

} } // namespace CORK::basis


namespace CORK {

  template <typename Sequence>
  struct value_type< basis::functions<Sequence> >
  : std::invoke_result< typename std::decay< Sequence >::type::value_type >
  //: std::result_of< typename std::decay< Sequence >::type::value_type >
  {
 //   typedef typename std::decay< Sequence >::type::value_type::result_type type ;
  } ;

  template <typename T, typename Sequence>
  struct value_type_for< T, basis::functions<Sequence> >
  : value_type< basis::functions<Sequence> >
  {} ;

} // namespace CORK

namespace std {
 
  template <typename T, typename S>
  struct invoke_result< std::function<T(S)> > {
    typedef T type ;
  } ;
}

#endif
