//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_barycentric_rational_hpp
#define cork_basis_barycentric_rational_hpp

#include <cork/utility/is_infinite.hpp>
#include <cork/concept/arithmetic.hpp>
#include <cork/concept/has_type_member.hpp>
#include <cork/basis/explicit_iterator.hpp>
#include <type_traits>
#include <cassert>
#define CORK_BARYCENTRIC_RATIONAL_WEIGHTED

namespace CORK { namespace basis {

  template <typename Weights, typename Nodes, typename I=typename std::decay< Weights >::type::size_type>
  class barycentric_rational
  {
    public:
      typedef typename std::decay< Weights >::type weights_type ;
      typedef typename std::decay< Nodes >::type   nodes_type ;
      typedef I                                    size_type ;

    public:
      explicit barycentric_rational( Weights weights, Nodes nodes )
      : weights_( weights )
      , nodes_( nodes )
      {
        assert( nodes_.size()==weights_.size() ) ;
      } // barycentric_rational

    public:
      I num_terms() const {
        return nodes_.size()+1 ;
      } // num_terms()

      nodes_type const&  nodes() const {
        return nodes_ ;
      } // nodes()

      typename std::add_lvalue_reference< Nodes >::type nodes() {
        return nodes_ ;
      } // nodes()

      weights_type const&  weights() const {
        return weights_ ;
      } // nodes()

      typename std::add_lvalue_reference< Weights >::type weights() {
        return weights_ ;
      } // nodes()

    private:
      typedef typename std::common_type< typename weights_type::value_type, typename nodes_type::value_type>::type inner_value_type ;

    public:
      typedef inner_value_type minimal_value_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename std::common_type< ValueType, inner_value_type >::type ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< T, inner_value_type > > ;

    public:
      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type_for<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;
        values(0) = 1.0 ;
        for (typename std::decay<FunctionValues>::type::size_type i=0; i<num_terms()-1; ++i) {
          if (is_infinite(nodes_(i))) {
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
            //values(i+1) = weights_(i) * arg ;
            values(i+1) = weights_(i) ;
#else
            //values(i+1) = arg ;
            values(i+1) = 1.0 ;
#endif
          } else {
            assert( !is_infinite(nodes_(i)) ) ;
            if (arg==nodes_(i)) {
              fill( values(glas2::range_from_end(1,0)), 0.0 ) ;
              values(i+1) = 1. ; /// weights_(i) ;
              return ;
            } else
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
              values(i+1) = weights_(i) / (arg - nodes_(i)) ;
#else
              values(i+1) = 1.0 / (arg - nodes_(i)) ;
#endif
          }
        }
#ifdef CORK_BARYCENTRIC_RATIONAL_WEIGHTED
        values(glas2::range_from_end(1,0)) /= sum(values(glas2::range_from_end(1,0))) ;
#else
        values(glas2::range_from_end(1,0)) /= sum(weights_ * values(glas2::range_from_end(1,0))) ;
#endif
      } // evaluate()

    public:
      template <typename ShiftValueType>
      using iterator = explicit_iterator< value_type_for<ShiftValueType> > ;

      template <typename ShiftValueType>
      decltype (auto) evaluate_iterator( ShiftValueType const& arg ) const {
        return explicit_iterator<ShiftValueType>( *this, arg ) ;
      } // evaluate_iterator

    private:
      Weights weights_ ;
      Nodes   nodes_ ;
  } ; // class barycentric_rational


  template <typename Weights, typename Nodes>
  barycentric_rational<Weights&,Nodes&> make_barycentric_rational_lvalue( Weights& weights, Nodes& nodes ) {
    return barycentric_rational<Weights&, Nodes&>( weights, nodes ) ;
  }

  template <typename Weights, typename Nodes>
  barycentric_rational<Weights&&,Nodes&&> make_barycentric_rational_rvalue( Weights&& weights, Nodes&& nodes ) {
    return barycentric_rational<Weights&&,Nodes&&>( weights, nodes ) ;
  }

  template <typename Weights, typename Nodes>
  barycentric_rational<Weights,Nodes> make_barycentric_rational( Weights const& weights, Nodes const& nodes ) {
    return barycentric_rational<Weights,Nodes>( weights, nodes ) ;
  }

} } // namespace CORK::basis

#endif
