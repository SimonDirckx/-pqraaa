//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_pair_hpp
#define cork_basis_pair_hpp

#include <cork/concept/has_type_member.hpp>
#include <cork/basis/explicit_iterator.hpp>
#include <cork/utility/value_type_for.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace basis {

  //
  // We assume that the first elements of Basis1 and Basis2 are the same
  //
  template <typename Basis1, typename Basis2>
  class pair
  {
    public:
      typedef typename std::decay<Basis1>::type basis_1_type ;
      typedef typename std::decay<Basis2>::type basis_2_type ;
      typedef typename std::common_type< typename basis_1_type::size_type
                                       , typename basis_2_type::size_type
                                       >::type                             size_type ;

    public:
      template <typename T>
      using value_type_for = typename std::conditional< basis_1_type::template has_value_type_for<T>::value
                                                      , typename basis_1_type::template value_type_for<T>
                                                      , typename basis_2_type::template value_type_for<T>
                                                      >::type ;

      template <typename T>
      using has_value_type_for = std::bool_constant< basis_1_type::template has_value_type_for<T>::value || basis_2_type::template has_value_type_for<T>::value > ;

    public:
      explicit pair( Basis1 basis1, Basis2 basis2 )
      : basis_1_( basis1 )
      , basis_2_( basis2 )
      {} // pair

    public:
      size_type num_terms() const {
        assert( basis_1_.num_terms()==basis_2_.num_terms() ) ;
        return basis_1_.num_terms() ;
      } // num_terms()

      basis_1_type const& basis_1() const { return basis_1_ ; }
      basis_2_type const& basis_2() const { return basis_2_ ; }

    private:
      template <typename ShiftValueType, typename FunctionValues>
      void evaluate_impl( ShiftValueType const& arg, FunctionValues values, std::true_type ) const {
        basis_1_.evaluate( arg, values ) ;
      }

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate_impl( ShiftValueType const& arg, FunctionValues values, std::false_type ) const {
        basis_2_.evaluate( arg, values ) ;
      }

    public:
      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type_for<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        evaluate_impl( arg, values, typename basis_1_type::template has_value_type_for<ShiftValueType>() ) ;
      } // evaluate

    private:
      Basis1 basis_1_ ;
      Basis2 basis_2_ ;
  } ; // class pair


  template <typename Basis1, typename Basis2>
  auto make_pair_lvalue( Basis1& basis1, Basis2& basis2 ) {
    return pair<Basis1&, Basis2&>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  auto make_pair_rvalue( Basis1&& basis1, Basis2&& basis2 ) {
    return pair<Basis1&&, Basis2&&>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  auto make_pair( Basis1 const& basis1, Basis2 const& basis2 ) {
    return pair<Basis1, Basis2>( basis1, basis2 ) ;
  }

} } // namespace CORK::basis

#endif
