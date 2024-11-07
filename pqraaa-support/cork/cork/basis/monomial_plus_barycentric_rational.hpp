//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_monomial_plus_barycentric_rational_hpp
#define cork_basis_monomial_plus_barycentric_rational_hpp

#include <cork/basis/monomial.hpp>
#include <glas2/vector.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace basis {

  template <typename BarycentricRational>
  class monomial_plus_barycentric_rational
  {
    public:
      typedef BarycentricRational barycentric_rational_type ;
      typedef typename barycentric_rational_type::size_type size_type ;

    public:
      explicit monomial_plus_barycentric_rational( size_type grade, barycentric_rational_type const& barycentric )
      : monomial_( grade )
      , barycentric_( barycentric )
      {}

    public:
      size_type num_terms() const {
        return monomial_.num_terms()+barycentric_.num_terms()-1 ;
      } // num_terms()

    public:
      template <typename ValueType>
      using value_type_for = typename barycentric_rational_type::template value_type_for<ValueType> ;

      template <typename ValueType>
      using has_value_type_for = typename barycentric_rational_type::template has_value_type_for<ValueType> ;

    public:
      // basis = [1, s, s^2, s^{d-1}, b_1,...,b_r] with b_i barycentric basis functions
      // The first grade+1 components of values are the values of the monomial part
      // The remaining components are the rational terms.
      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        assert( values.size() == num_terms() ) ;

        glas2::vector< ShiftValueType > mono( monomial_.num_terms() ) ;
        monomial_.evaluate( arg, mono ) ;

        glas2::vector< ShiftValueType > bary( barycentric_.num_terms() ) ;
        barycentric_.evaluate( arg, bary ) ;

        values( glas2::range(0,monomial_.grade()+1) ) = mono ;
        values( glas2::range_from_end(monomial_.grade()+1,0) ) = bary(glas2::range_from_end(1,0)) ;
      } // evaluate()

    public:
      auto const& monomial() const { return monomial_ ; }
      auto const& barycentric() const { return barycentric_ ; }

    private:
      basis::monomial<size_type>   monomial_ ;
      barycentric_rational_type    barycentric_ ;
  } ; // class monomial_plus_barycentric_rational


} } // namespace CORK::basis

#endif
