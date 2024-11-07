//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_accumulated_partial_fractions_real_hpp
#define cork_basis_accumulated_partial_fractions_real_hpp

#include <cork/concept/has_type_member.hpp>
#include <cork/basis/barycentric_rational_real.hpp>
#include <type_traits>
#include <cassert>
#include <complex>

namespace CORK { namespace basis {

  template <typename Weights, typename Poles, typename I=typename std::decay< Poles >::type::size_type>
  class accumulated_partial_fractions_real
  {
    public:
      typedef typename std::decay< Weights >::type weights_type ;
      typedef typename std::decay< Poles >::type   poles_type ;
      typedef I                                    size_type ;

    public:
      explicit accumulated_partial_fractions_real( Weights weights, Poles poles )
      : poles_( poles )
      , weights_( weights )
      , num_real_( 0 )
      {
        // Real poles come first
#ifndef NDEBUG
        bool already_complex = false ;
#endif
        assert( poles_.size()==weights_.size() ) ;
        for (size_type i=0; i<poles_.size(); ++i) {
          if (poles_(i).imag()!=0.0) {
            assert( poles_(i).real()==poles_(i+1).real() ) ;
            assert( std::real(weights_(i))==std::real(weights_(i+1)) ) ;
            assert( poles_(i).imag()==-poles_(i+1).imag() ) ;
            assert( std::imag(weights(i))==-std::imag(weights(i+1)) ) ;
#ifndef NDEBUG
            already_complex = true ;
#endif
            ++i ;
          } else {
            assert(!already_complex) ;
            assert( std::imag(weights_(i))==0.0 ) ;
            ++num_real_ ;
          }
        }
      } // accumulated_partial_fractions_real

    public:
      I num_terms() const {
        return poles_.size()+1 ;
      } // num_terms()

      I num_real() const { return num_real_ ; }

      poles_type const&  poles() const { return poles_ ; } // poles()
      weights_type const&  weights() const { return weights_ ; } // weights()

      typename std::add_lvalue_reference< Poles >::type poles() { return poles_ ; } // poles()
      typename std::add_lvalue_reference< Weights >::type weights() { return weights_ ; } // weights()

    public:
      typedef typename std::common_type< typename weights_type::value_type, typename poles_type::value_type >::type data_value_type ;

    private:
      typedef decltype(std::abs(data_value_type())) real_data_value_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename std::conditional< std::is_arithmetic<ValueType>::value
                                                      , typename std::common_type< ValueType, real_data_value_type >::type
                                                      , typename std::common_type< ValueType, data_value_type >::type
                                                      >::type ;

      template <typename T>
      using has_value_type_for = std::bool_constant< has_type_member< std::common_type< T, real_data_value_type > >::value || has_type_member< std::common_type< T, data_value_type > >::value > ;

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type_for<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;
        values(0) = 1.0 ;
        typename std::decay<decltype( weights_(0) )>::type previous_weight = 1.0 ;
        data_value_type previous_value1 = 1.0 ;
        data_value_type previous_value2 = 1.0 ;

        // Real poles
        for (typename std::decay<FunctionValues>::type::size_type i=0; i<num_real_; ++i) {
          assert( arg!=poles_(i) ) ;
          assert( poles_(i).imag()==0. ) ;
          assert( std::imag(weights_(i))==0. ) ;
          barycentric_rational_real_detail::convert_to(values(i+1)) = std::real(weights_(i)) * previous_value1 / ( std::real(previous_weight) * (arg - poles_(i).real()) ) ;
          previous_weight = weights_(i) ;
          previous_value1 = previous_value2 = values(i+1) ;
        }

        // Complex poles
        for (typename std::decay<FunctionValues>::type::size_type i=num_real_; i<num_terms()-1; i+=2) {
          assert( std::conj(poles_(i))==poles_(i+1) ) ;
          assert( std::conj(weights_(i))==weights_(i+1) ) ;

          auto new_value1 = previous_value1 * weights_(i)/ ( previous_weight * (arg - poles_(i)) ) ;
          auto new_value2 = previous_value2 * weights_(i+1)/ ( std::conj(previous_weight) * (arg - poles_(i+1)) ) ;
          barycentric_rational_real_detail::convert_to(values(i+1)) = new_value1 + new_value2 ;
          barycentric_rational_real_detail::convert_to(values(i+2)) = ( new_value1 - new_value2 ) * data_value_type(0.,1.) ;

          previous_weight = weights_(i) ;
          previous_value1 = new_value1 ;
          previous_value2 = new_value2 ;
        }
      } // evaluate

/*    public:
      template <typename ShiftValueType>
      using iterator = explicit_iterator< value_type<ShiftValueType> > ;

      template <typename ShiftValueType>
      decltype (auto) evaluate_iterator( ShiftValueType const& arg ) const {
        return explicit_iterator<ShiftValueType>( *this, arg ) ;
      } // evaluate_iterator
*/
    private:
      Poles     poles_ ;
      Weights   weights_ ;
      size_type num_real_ ;
  } ; // class accumulated_partial_fractions_real


  template <typename Poles, typename Weights>
  auto make_accumulated_partial_fractions_real_lvalue( Weights& weights, Poles& poles ) {
    return accumulated_partial_fractions_real<Weights&,Poles&>( weights, poles ) ;
  }

  template <typename Poles, typename Weights>
  decltype (auto) make_accumulated_partial_fractions_real_rvalue( Weights&& weights, Poles&& poles ) {
    return accumulated_partial_fractions_real<Weights&&, Poles&&>( weights, poles ) ;
  }

  template <typename Poles, typename Weights>
  auto make_accumulated_partial_fractions_real( Weights const& weights, Poles const& poles ) {
    return accumulated_partial_fractions_real<Weights, Poles>( weights, poles ) ;
  }

} } // namespace CORK::basis

#endif
