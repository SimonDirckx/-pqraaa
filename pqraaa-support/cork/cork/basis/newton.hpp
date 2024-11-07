//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_newton_hpp
#define cork_basis_newton_hpp

#include <cork/basis/iterator.hpp>
#include <type_traits>
#include <cassert>

namespace CORK { namespace basis {

  template <typename Points, typename I=typename std::decay< Points >::type::size_type>
  class newton
  {
    public:
      typedef typename std::decay< Points >::type  points_type ;
      typedef I                                    size_type ;

    public:
      explicit newton( Points points )
      : points_( points )
      {} // newton

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

    public:
      template <typename ShiftValueType>
      using value_type = typename std::common_type< typename points_type::value_type, ShiftValueType >::type ;

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_assert( std::is_convertible< value_type<ShiftValueType>, typename std::decay<FunctionValues>::type::value_type >::value, "" ) ;
        assert( values.size() == num_terms() ) ;
        values(0) = 1.0 ;
        for (typename std::decay<FunctionValues>::type::size_type i=1; i<values.size(); ++i) {
          values(i) = (arg - points_(i-1)) * values(i-1) ;
        }
      } // evaluate

    public:
      template <typename ShiftValueType>
      struct functor_type {
        functor_type( ShiftValueType const& arg, points_type const& points )
        : arg_( arg )
        , points_( points )
        , index_( 0 )
        {}

        template <typename ValueType>
        decltype (auto) operator() ( ValueType const& v ) {
          ++index_ ;
          return v * (arg_ - points_(index_-1)) ;
        }

        ShiftValueType     arg_ ;
        points_type const& points_ ;
        size_type          index_ ;
      } ;

      template <typename ShiftValueType>
      using iterator = CORK::basis::iterator< value_type<ShiftValueType>, functor_type<ShiftValueType> > ;

      template <typename ShiftValueType>
      decltype (auto) evaluate_iterator( ShiftValueType const& arg ) const {
        return iterator<ShiftValueType>( functor_type<ShiftValueType>( arg, points_ ) ) ;
      } // evaluate_iterator

    private:
      Points points_ ;

  } ; // class newton


  template <typename Points>
  newton<Points&> make_newton_lvalue( Points& points ) {
    return newton<Points&>( points ) ;
  }

  template <typename Points>
  newton<Points&&> make_newton_rvalue( Points&& points ) {
    return newton<Points&&>( points ) ;
  }

  template <typename Points>
  newton<Points> make_newton( Points const& points ) {
    return newton<Points>( points ) ;
  }

} } // namespace CORK::basis

#endif
