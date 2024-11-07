//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_union_of_functions_hpp
#define cork_basis_union_of_functions_hpp

#include <cork/concept/has_type_member.hpp>
#include <cork/basis/concept.hpp>
#include <cork/basis/explicit_iterator.hpp>
#include <cork/utility/value_type_for.hpp>
#include <cork/utility/ref.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace basis {

  //
  // We assume that the first elements of Basis1 and Basis2 are the same
  //
#ifdef CORK_HAS_CONCEPTS
  template <Concept::Basis Basis1, Concept::Basis Basis2>
#else
  template <typename Basis1, typename Basis2>
#endif
  class union_of_functions
  {
    public:
      typedef typename CORK::deref_type<Basis1>::type basis_1_type ;
      typedef typename CORK::deref_type<Basis2>::type basis_2_type ;
      typedef typename std::common_type< typename basis_1_type::size_type
                                       , typename basis_2_type::size_type
                                       >::type                             size_type ;

    public:
      template <typename T>
      using value_type_for = typename std::common_type< typename basis_1_type::template value_type_for<T>, typename basis_2_type::template value_type_for<T> >::type ;

      template <typename T>
      using has_value_type_for = std::bool_constant< basis_1_type::template has_value_type_for<T>::value && basis_2_type::template has_value_type_for<T>::value > ;

    public:
      explicit union_of_functions( Basis1 const& basis1, Basis2 const& basis2 )
      : basis_1_( basis1 )
      , basis_2_( basis2 )
      {} // union_of_functions

    public:
      size_type num_terms() const {
        return CORK::deref(basis_1_).num_terms() + CORK::deref(basis_2_).num_terms() ;
      } // num_terms()

      basis_1_type const& basis_1() const { return CORK::deref(basis_1_) ; }
      basis_2_type const& basis_2() const { return CORK::deref(basis_2_) ; }

      Basis1 const& basis_1_storage() const { return basis_1_ ; }
      Basis2 const& basis_2_storage() const { return basis_2_ ; }

    public:
      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type_for<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( (long int)(values.size()) == (long int)(num_terms()) ) ;

        CORK::deref(basis_2_).evaluate( arg, values( glas2::range_from_end( CORK::deref(basis_1_).num_terms(), 0 ) ) ) ;
        CORK::deref(basis_1_).evaluate( arg, values( glas2::range( 0, CORK::deref(basis_1_).num_terms() ) ) ) ;
      } // evaluate

    private:
      Basis1 basis_1_ ;
      Basis2 basis_2_ ;
  } ; // class union_of_functions


  template <typename Basis1, typename Basis2>
  union_of_functions<Basis1&, Basis2&> make_union_of_functions_lvalue( Basis1& basis1, Basis2& basis2 ) {
    return union_of_functions<Basis1&, Basis2&>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  union_of_functions<Basis1&&, Basis2&&> make_union_of_functions_rvalue( Basis1&& basis1, Basis2&& basis2 ) {
    return union_of_functions<Basis1&&, Basis2&&>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  union_of_functions<Basis1, Basis2> make_union_of_functions( Basis1 const& basis1, Basis2 const& basis2 ) {
    return union_of_functions<Basis1, Basis2>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  std::ostream& operator<<( std::ostream& s, union_of_functions<Basis1,Basis2> const& b ) {
    s << "union_of_bases:{" << b.basis_1() << "," << b.basis_2() << "}" ;
    return s ;
  } // operator<<

} } // namespace CORK::basis

#endif
