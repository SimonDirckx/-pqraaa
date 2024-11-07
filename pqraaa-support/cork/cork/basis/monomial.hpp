//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_monomial_hpp
#define cork_basis_monomial_hpp

#include <cork/concept/arithmetic.hpp>
#include <cork/basis/iterator.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace basis {

  template <typename I=int>
  class monomial
  {
    public:
      typedef I size_type ;
      static_assert( std::is_integral<I>::value, "basis::monomial<I>: I should be integral type" ) ;

    public:
      explicit monomial( size_type grade )
      : grade_( grade )
      {
        assert( grade_ >= 0 ) ;
      } // monomial

    public:
      size_type num_terms() const {
        return grade_+1 ;
      } // num_terms

      size_type grade() const {
        return grade_ ;
      } // grade

    public:
      template <typename ValueType, typename FunctionValues>
      void evaluate( ValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< ValueType, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;
        values(0) = 1.0 ;
        for (typename std::decay<FunctionValues>::type::size_type i=1; i<values.size(); ++i) {
          values(i) = arg * values(i-1) ;
        }
      } // evaluate

    public:
      template <typename ValueType>
      using value_type_for = ValueType ;

      template <typename T>
      using has_value_type_for = CORK::is_arithmetic<T> ;

      template <typename T>
      using has_compatible_argument_type = std::true_type ;

      template <typename T>
      using has_compatible_result_type = std::true_type ;

    public:
      template <typename ValueType>
      struct functor_type {
        functor_type( ValueType const& arg )
        : arg_( arg )
        {}

        ValueType operator() ( ValueType const& v ) const {
          return v * arg_ ;
        }

        ValueType arg_ ;
      } ;

      template <typename ValueType>
      using iterator = CORK::basis::iterator<ValueType, functor_type<ValueType> > ;

      template <typename ValueType>
      decltype(auto) evaluate_iterator( ValueType const& arg ) const {
        return iterator<ValueType>( functor_type<ValueType>( arg ) ) ;
      } // evaluate_iterator


    private:
      size_type grade_ ;

  } ; // class monomial

  template <typename I>
  std::ostream& operator<<( std::ostream& s, monomial<I> const& b ) {
    s << "monomial:{degree=" << b.grade() << "}" ;
    return s ;
  } // operator<<

} } // namespace CORK::basis

#endif
