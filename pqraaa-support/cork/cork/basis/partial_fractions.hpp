//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_partial_fractions_hpp
#define cork_basis_partial_fractions_hpp

#include <cork/basis/explicit_iterator.hpp>
#include <cork/concept/has_type_member.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace basis {

  template <typename Poles, typename Weights, typename I=typename std::decay< Poles >::type::size_type>
  class partial_fractions
  {
    public:
      typedef typename std::decay< Poles >::type   poles_type ;
      typedef typename std::decay< Weights >::type weights_type ;
      typedef I                                    size_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename std::common_type< ValueType, typename poles_type::value_type, typename weights_type::value_type >::type ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< T, typename poles_type::value_type, typename weights_type::value_type > > ;

    public:
      explicit partial_fractions( Poles poles, Weights weights )
      : poles_( poles )
      , weights_( weights )
      {
        assert( poles_.size()==weights_.size() ) ;
      } // partial_fractions

    public:
      I num_terms() const {
        return poles_.size()+1 ;
      } // num_terms()

      poles_type const&  poles() const { return poles_ ; } // poles()
      weights_type const&  weights() const { return weights_ ; } // weights()

      typename std::add_lvalue_reference< Poles >::type poles() { return poles_ ; } // poles()
      typename std::add_lvalue_reference< Weights >::type weights() { return weights_ ; } // weights()

    private:
      typedef typename poles_type::value_type inner_value_type ;

    public:
      template <typename ShiftValueType>
      using value_type = typename std::common_type< inner_value_type, ShiftValueType >::type ;

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;
        values(0) = 1.0 ;
        for (typename std::decay<FunctionValues>::type::size_type i=1; i<num_terms(); ++i) {
          assert( arg!=poles_(i-1) ) ;
          values(i) = weights_(i-1) / (arg - poles_(i-1)) ;
        }
      } // evaluate

    public:
      template <typename ShiftValueType>
      using iterator = explicit_iterator< value_type<ShiftValueType> > ;

      template <typename ShiftValueType>
      decltype (auto) evaluate_iterator( ShiftValueType const& arg ) const {
        return explicit_iterator<ShiftValueType>( *this, arg ) ;
      } // evaluate_iterator

    private:
      Poles   poles_ ;
      Weights weights_ ;
  } ; // class partial_fractions


  template <typename Poles, typename Weights>
  partial_fractions<Poles&,Weights&> make_partial_fractions_lvalue( Poles& poles, Weights& weights ) {
    return partial_fractions<Poles&,Weights&>( poles, weights ) ;
  }

  template <typename Poles, typename Weights>
  decltype (auto) make_partial_fractions_rvalue( Poles&& poles, Weights&& weights ) {
    return partial_fractions<Poles&&, Weights&&>( poles, weights ) ;
  }

  template <typename Poles, typename Weights>
  partial_fractions<Poles, Weights> make_partial_fractions( Poles const& poles, Weights const& weights ) {
    return partial_fractions<Poles, Weights>( poles, weights ) ;
  }

} } // namespace CORK::basis

#endif
