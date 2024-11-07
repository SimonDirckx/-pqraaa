//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_inverse_newton_hpp
#define cork_basis_inverse_newton_hpp

#include <cork/basis/iterator.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace basis {

  template <typename Points, typename Weights, typename I=typename std::decay< Points >::type::size_type>
  class inverse_newton
  {
    public:
      typedef typename std::decay< Points >::type  points_type ;
      typedef typename std::decay< Weights >::type weights_type ;
      typedef I                                    size_type ;

    public:
      explicit inverse_newton( Points points, Weights weights )
      : points_( points )
      , weights_( weights )
      {} // inverse_newton

    public:
      I num_terms() const {
        return points_.size()+1 ;
      } // num_terms

      points_type const&  points() const {
        return points_ ;
      } // point

      typename std::add_lvalue_reference< Points >::type points() {
        return points_ ;
      } // points

      weights_type const&  weights() const {
        return weights_ ;
      } // weights

      typename std::add_lvalue_reference< Weights >::type weights() {
        return weights_ ;
      } // weights

    public:
      template <typename ShiftValueType>
      using value_type = typename std::common_type< typename std::common_type<typename points_type::value_type, typename weights_type::value_type>::type, ShiftValueType >::type ;

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        typedef typename FunctionValues::value_type value_type ;

        assert( values.size() == num_terms() ) ;
        values(0) = 1.0 ;
        for (typename std::decay<FunctionValues>::type::size_type i=1; i<values.size(); ++i) {
          assert( glas2::imag(weights(i-1))==0.0 ) ;
          if (glas2::imag(points_(i-1))==0.0) {
            values(i) = weights_(i) / (arg - points_(i-1)) * values(i-1) ;
          } else {
            value_type val1 = weights_(i) / (arg - points_(i-1)) * values(i-1) ;
            value_type val2 = weights_(i) / (arg - points_(i)) * values(i-1) ;
            values(i) = 0.5*(val1+val2) ;
            values(i+1) = value_type(0.,-0.5)*(val1+val2) ;
            ++i ;
          }
        }
      } // evaluate

    private:
      Points  points_ ;
      Weights weights_ ;

  } ; // class inverse_newton


  template <typename Points>
  inverse_newton<Points&> make_inverse_newton_lvalue( Points& points ) {
    return inverse_newton<Points&>( points ) ;
  }

  template <typename Points>
  inverse_newton<Points&&> make_inverse_newton_rvalue( Points&& points ) {
    return inverse_newton<Points&&>( points ) ;
  }

  template <typename Points>
  inverse_newton<Points> make_inverse_newton( Points const& points ) {
    return inverse_newton<Points>( points ) ;
  }

} } // namespace CORK::basis

#endif
