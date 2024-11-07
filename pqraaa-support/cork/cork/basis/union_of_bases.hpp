//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_union_of_bases_hpp
#define cork_basis_union_of_bases_hpp

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
  class union_of_bases
  {
    public:
      typedef typename std::decay<Basis1>::type basis_1_type ;
      typedef typename std::decay<Basis2>::type basis_2_type ;
      typedef typename std::common_type< typename basis_1_type::size_type
                                       , typename basis_2_type::size_type
                                       >::type                             size_type ;

    public:
      explicit union_of_bases( Basis1 basis1, Basis2 basis2 )
      : basis_1_( basis1 )
      , basis_2_( basis2 )
      {} // union_of_bases

    public:
      size_type num_terms() const {
        return basis_1_.num_terms() + basis_2_.num_terms() - 1 ;
      } // num_terms()

      basis_1_type const& basis_1() const { return basis_1_ ; }
      basis_2_type const& basis_2() const { return basis_2_ ; }

      basis_1_type& basis_1() { return basis_1_ ; }
      basis_2_type& basis_2() { return basis_2_ ; }

    public:
      template <typename ShiftValueType>
      using value_type_for = typename std::common_type< typename basis_1_type::template value_type_for< ShiftValueType >
                                                      , typename basis_2_type::template value_type_for< ShiftValueType >
                                                      >::type ;

      template <typename ShiftValueType>
      using has_value_type_for = std::bool_constant< basis_1_type::template has_value_type_for<ShiftValueType>::value && basis_2_type::template has_value_type_for<ShiftValueType>::value > ;

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type_for<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;

        basis_2_.evaluate( arg, values( glas2::range_from_end( basis_1_.num_terms()-1, 0 ) ) ) ;
#ifndef NDEBUG
        value_type_for<ShiftValueType> val0 = values(basis_1_.num_terms()-1) ;
#endif
        basis_1_.evaluate( arg, values( glas2::range( 0, basis_1_.num_terms() ) ) ) ;
#ifndef NDEBUG
        assert( values(0) == val0 ) ;
#endif
      } // evaluate

    public:
      template <typename ShiftValueType>
      using iterator = explicit_iterator< value_type_for<ShiftValueType> > ;

      template <typename ShiftValueType>
      decltype (auto) evaluate_iterator( ShiftValueType const& arg ) const {
        return explicit_iterator<ShiftValueType>( *this, arg ) ;
      } // evaluate_iterator

    private:
      Basis1 basis_1_ ;
      Basis2 basis_2_ ;
  } ; // class union_of_bases


  template <typename Basis1, typename Basis2>
  union_of_bases<Basis1&, Basis2&> make_union_of_bases_lvalue( Basis1& basis1, Basis2& basis2 ) {
    return union_of_bases<Basis1&, Basis2&>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  union_of_bases<Basis1&&, Basis2&&> make_union_of_bases_rvalue( Basis1&& basis1, Basis2&& basis2 ) {
    return union_of_bases<Basis1&&, Basis2&&>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  union_of_bases<Basis1, Basis2> make_union_of_bases( Basis1 const& basis1, Basis2 const& basis2 ) {
    return union_of_bases<Basis1, Basis2>( basis1, basis2 ) ;
  }

  template <typename Basis1, typename Basis2>
  std::ostream& operator<<( std::ostream& s, union_of_bases<Basis1,Basis2> const& b ) {
    s << "union_of_bases:{" << b.basis_1() << "," << b.basis_2() << "}" ;
    return s ;
  } // operator<<

} } // namespace CORK::basis

#endif
