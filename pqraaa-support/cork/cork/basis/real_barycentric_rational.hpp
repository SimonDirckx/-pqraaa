//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_real_barycentric_rational_hpp
#define cork_basis_real_barycentric_rational_hpp

#include <cork/basis/explicit_iterator.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace basis {

  //
  // Poles must be complex valued
  //
  template <typename Weights, typename Poles, typename I=typename std::decay< Weights >::type::size_type>
  class real_barycentric_rational
  {
    public:
      typedef typename std::decay< Weights >::type weights_type ;
      typedef typename std::decay< Poles >::type   poles_type ;
      typedef I                                    size_type ;

    public:
      explicit real_barycentric_rational( Weights weights, Poles poles )
      : weights_( weights.size() )
      , poles_real_( poles.size() )
      , poles_imag_( poles.size() )
      {
        poles_real_ = glas2::real(poles) ;
        poles_imag_ = glas2::imag(poles) ;
        weights_ = real(weights) ;

        for (size_type i=0; i<weights.size(); ++i) {
          if (real(poles(i))!=0.0) {
            assert( poles(i)==std::conj(poles(i+1)) ) ;
            assert( weights(i)==std::conj(weights(i+1)) ) ;
            weights_(i+1) = imag(weights(i)) ;
            ++i ;
        }
      } // real_barycentric_rational

    public:
      I num_terms() const {
        return weights_.size()+1 ;
      } // num_terms()

      poles_real_type const&  poles_real() const {
        return poles_real_ ;
      } // poles_real()

      poles_imag_type const&  poles_imag() const {
        return poles_imag_ ;
      } // poles_real()

      weights_type const&  weights() const {
        return weights_ ;
      } // poles()

    private:
      typedef typename std::common_type< typename weights_type::value_type, typename poles_type::value_type::>::type inner_value_type ;

    public:
      template <typename ShiftValueType>
      using value_type = typename std::common_type< inner_value_type, ShiftValueType >::type ;

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;
        values(0) = 1.0 ;
        for (typename std::decay<FunctionValues>::type::size_type i=0; i<num_terms()-1; ++i) {
          assert( arg!=poles_(i) ) ;
          if (imag(poles_(i))==0)
            values(i+1) = 1.0 / (arg - real(poles_(i))) ;
          else {
            assert(i<num_terms()-2) ;
            assert( imag(poles_(i))==-imag(poles_(i+1)) ) ;
            values(i+1) = (arg - real(poles_(i))) / (arg*arg - 2*arg* real(poles_(i)) + std::norm(poles_(i))) ;
            values(i+2) = imag(poles_(i)) / (arg*arg - 2*arg* real(poles_(i)) + std::norm(poles_(i))) ;
          }
        }
        values(glas2::range_from_end(1,0)) /= sum(weights_*values(glas2::range_from_end(1,0))) ;
      } // evaluate

    public:
      template <typename ShiftValueType>
      using iterator = explicit_iterator< value_type<ShiftValueType> > ;

      template <typename ShiftValueType>
      decltype (auto) evaluate_iterator( ShiftValueType const& arg ) const {
        return explicit_iterator<ShiftValueType>( *this, arg ) ;
      } // evaluate_iterator

    private:
      Weights weights_ ;
      Poles   poles_ ;
  } ; // class real_barycentric_rational


  template <typename Weights, typename Poles>
  real_barycentric_rational<Weights&,Poles&> make_real_barycentric_rational_lvalue( Weights& weights, Poles& poles ) {
    return real_barycentric_rational<Weights&, Poles&>( weights, poles ) ;
  }

  template <typename Weights, typename Poles>
  real_barycentric_rational<Weights&&,Poles&&> make_real_barycentric_rational_rvalue( Weights&& weights, Poles&& poles ) {
    return real_barycentric_rational<Weights&&,Poles&&>( weights, poles ) ;
  }

  template <typename Weights, typename Poles>
  real_barycentric_rational<Weights,Poles> make_real_barycentric_rational( Weights const& weights, Poles const& poles ) {
    return real_barycentric_rational<Weights,Poles>( weights, poles ) ;
  }

} } // namespace CORK::basis

#endif
