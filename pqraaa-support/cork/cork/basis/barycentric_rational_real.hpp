//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_barycentric_rational_real_hpp
#define cork_basis_barycentric_rational_real_hpp

#include <cork/concept/has_type_member.hpp>
#include <cork/utility/is_infinite.hpp>
#include <cork/basis/explicit_iterator.hpp>
#include <type_traits>
#include <cassert>
#include <cmath>
#include <limits>

namespace CORK { namespace basis {

  namespace barycentric_rational_real_detail {
    template <typename T, typename C, typename EnableIf=void>
    struct select {
    } ;

    template <typename T>
    struct select< T, T, void > {
      static void assign( T& to, T const& from ) { to = from ; }

      static void plus_assign( T& to, T const& from ) { to += from ; }

      static void minus_assign( T& to, T const& from ) { to -= from ; }
    } ;

    template <typename R>
    struct select< std::complex<R>, R, typename std::enable_if<std::is_arithmetic<R>::value >::type > {
      static void assign( std::complex<R>& to, R const& from ) { to = from ; }

      static void plus_assign( std::complex<R>& to, R const& from ) { to += from ; }

      static void minus_assign( std::complex<R>& to, R const& from ) { to -= from ; }
    } ;

    template <typename R>
    struct select< R, std::complex<R>, typename std::enable_if<std::is_arithmetic<R>::value >::type > {
      static void assign( R& to, std::complex<R> const& from ) { assert( from.imag()==0.0) ; to = from.real() ; }

      static void plus_assign( R& to, std::complex<R> const& from ) { assert( from.imag()==0.0) ; to += from.real() ; }

      static void minus_assign( R& to, std::complex<R> const& from ) { assert( from.imag()==0.0) ; to -= from.real() ; }
    } ;

    template <typename T>
    class convert_to_type {
      public:
        convert_to_type( T& value )
        : value_( value )
        {}

        template <typename C>
        void operator=( C const& v ) {
          select<T,C>::assign( value_, v ) ;
        }

        template <typename C>
        void operator+=( C const& v ) {
          select<T,C>::plus_assign( value_, v ) ;
        }

        template <typename C>
        void operator-=( C const& v ) {
          select<T,C>::minus_assign( value_, v ) ;
        }

      private:
        T& value_ ;
    } ;

    template <typename T>
    convert_to_type<T> convert_to( T& r ) {
      return convert_to_type<T>( r ) ;
    }
  } // namespace barycentric_rational_real_detail

  template <typename Weights, typename Nodes, typename I=typename std::decay< Weights >::type::size_type>
  class barycentric_rational_real
  {
    public:
      typedef typename std::decay< Weights >::type weights_type ;
      typedef typename std::decay< Nodes >::type   nodes_type ;
      typedef I                                    size_type ;
      typedef typename std::common_type< typename weights_type::value_type, typename nodes_type::value_type >::type data_value_type ;

    public:
      explicit barycentric_rational_real( Weights weights, Nodes nodes )
      : weights_( weights )
      , nodes_( nodes )
      {
        // Real nodes come first
#ifndef NDEBUG
        bool already_complex = false ;
#endif
        assert( nodes_.size()==weights_.size() ) ;
        for (size_type i=0; i<nodes_.size(); ++i) {
          if (nodes_(i).imag()!=0.0) {
            assert( nodes_(i).real()==nodes_(i+1).real() ) ;
            assert( weights_(i).real()==weights_(i+1).real() ) ;
            assert( nodes_(i).imag()==-nodes_(i+1).imag() ) ;
            assert( weights(i).imag()==-weights(i+1).imag() ) ;
#ifndef NDEBUG
            already_complex = true ;
#endif
            ++i ;
          } else {
            assert(!already_complex) ;
            assert( weights_(i).imag()==0.0 ) ;
          }
        }
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
      typedef decltype(std::abs(inner_value_type())) real_inner_value_type ;

    public:
      typedef real_inner_value_type minimal_value_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename std::conditional< std::is_arithmetic<ValueType>::value
                                                      , typename std::common_type< ValueType, real_inner_value_type >::type
                                                      , typename std::common_type< ValueType, inner_value_type >::type
                                                      >::type ;

      template <typename T>
      using has_value_type_for = std::bool_constant< has_type_member< std::common_type< T, real_inner_value_type > >::value || has_type_member< std::common_type< T, inner_value_type > >::value > ;

    public:
      // If support point is finite, the basis function is proportional to weight(i)/(z-nodes(i))
      // If the support point is infinite, the basis function is proportional to weight(i)*z
      // Since we have to put z in M-z*N, we use the basis function weight(i) instead.
      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible<ShiftValueType, typename FunctionValues::value_type>::value
                     , "CORK::basis4CORK::barycentric_rational_real: ShiftValueType must be convertible to the value_type of FunctionValues"
                     ) ;
        assert( values.size() == num_terms() ) ;
        typedef std::complex< decltype(std::abs(typename FunctionValues::value_type())) > complex_value_type ;

        values(0) = 1.0 ;
        typename FunctionValues::value_type sum = 0.0 ;
        for (typename std::decay<FunctionValues>::type::size_type i=0; i<num_terms()-1; ++i) {
         if (nodes_(i).imag()==0.) {
            assert( weights_(i).imag()==0. ) ;
            if (arg==nodes_(i)) {
            //if (std::abs(arg-nodes_(i))<std::sqrt(std::numeric_limits<decltype(std::abs(arg))>::epsilon())) {
              fill( values(glas2::range_from_end(1,0)), 0.0 ) ;
              values(i+1) = 1. ; /// weights_(i) ;
              return ;
            } else {
              if (is_infinite(nodes_(i))) {
                barycentric_rational_real_detail::convert_to(values(i+1)) = weights_(i).real() ;
              } else {
                barycentric_rational_real_detail::convert_to(values(i+1)) = weights_(i).real() /(arg - nodes_(i).real()) ;
              }
              sum += values(i+1) ;
            }
          } else {
            assert( std::abs( weights_(i) - std::conj(weights_(i+1)))==0. ) ;
            assert( std::abs( nodes_(i) - std::conj(nodes_(i+1)))==0. ) ;
            if (arg==nodes_(i)) {
            //if (std::abs(arg-nodes_(i))<std::sqrt(std::numeric_limits<decltype(std::abs(arg))>::epsilon())) {
              fill( values(glas2::range_from_end(1,0)), 0.0 ) ;
              values(i+1) = 1. ;
              barycentric_rational_real_detail::convert_to(values(i+2)) = complex_value_type(0.,1.0) ;
              return ;
            } else if (arg==conj(nodes_(i))) {
            //} else if (std::abs(arg-conj(nodes_(i)))<std::sqrt(std::numeric_limits<decltype(std::abs(arg))>::epsilon())) {
              fill( values(glas2::range_from_end(1,0)), 0.0 ) ;
              values(i+1) = 1. ;
              barycentric_rational_real_detail::convert_to(values(i+2)) = complex_value_type(0.,-1.0) ;
              return ;
            } else {
              auto value_i = weights_(i) / (arg - nodes_(i)) ;
              auto value_i1 = weights_(i+1) / (arg - nodes_(i+1)) ;
              barycentric_rational_real_detail::convert_to(values(i+1)) = value_i + value_i1 ;
              sum += values(i+1) ;
              ++i ;
              barycentric_rational_real_detail::convert_to(values(i+1)) = complex_value_type(0.,1.0) * ( value_i - value_i1 ) ;
            }
          }
        }
        values(glas2::range_from_end(1,0)) /= sum ;
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
  } ; // class barycentric_rational_real


  template <typename Weights, typename Nodes>
  barycentric_rational_real<Weights&,Nodes&> make_barycentric_rational_real_lvalue( Weights& weights, Nodes& nodes ) {
    return barycentric_rational_real<Weights&, Nodes&>( weights, nodes ) ;
  }

  template <typename Weights, typename Nodes>
  barycentric_rational_real<Weights&&,Nodes&&> make_barycentric_rational_real_rvalue( Weights&& weights, Nodes&& nodes ) {
    return barycentric_rational_real<Weights&&,Nodes&&>( weights, nodes ) ;
  }

  template <typename Weights, typename Nodes>
  barycentric_rational_real<Weights,Nodes> make_barycentric_rational_real( Weights const& weights, Nodes const& nodes ) {
    return barycentric_rational_real<Weights,Nodes>( weights, nodes ) ;
  }

} } // namespace CORK::basis

#endif
