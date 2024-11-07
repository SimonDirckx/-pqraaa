//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_partial_fractions_hpp
#define cork_basis4cork_partial_fractions_hpp

#include <cork/basis/partial_fractions.hpp>
#include <cork/basis4cork/partial_fractions_handle.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <cassert>
#include <string>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename Poles, typename Weights, typename I>
  class basis4CORK< basis::partial_fractions<Poles,Weights,I> >
  {
    public:
      typedef typename std::decay<Poles>::type   poles_type ;
      typedef typename std::decay<Weights>::type weights_type ;
      typedef I                                  size_type ;

    public:
      template <typename ValueType>
      using value_type_for = typename std::common_type< ValueType, typename std::common_type< typename poles_type::value_type, typename weights_type::value_type>::type >::type ;

    public:
      explicit basis4CORK( basis::partial_fractions<Poles,Weights,I> const& basis )
      : poles_( basis.poles() )
      , weights_( basis.weights() )
      {}

    public:
      // Basis4CORK
      int num_terms() const { return poles_.size()+1 ; }
      int size() const { return num_terms() ; }

      typename std::add_lvalue_reference< typename std::add_const<Poles>::type >::type poles() const& { return poles_ ; }
      typename std::add_lvalue_reference< typename std::add_const<Weights>::type >::type weights() const& { return weights_ ; }

    public:
      // Basis4CORK
      template <typename UM, typename ZM>
      void multiply_N( UM const& U, ZM Z ) const {
        assert( U.num_rows()==size() ) ;
        assert( Z.num_rows()==size()-1 ) ;

       // fill( Z(0,glas2::all()), 0.0 ) ;
        Z = U( glas2::range(1,U.num_rows()), glas2::all() ) ;
      } // multiply_N()

      template <typename UM, typename ZM>
      void multiply_M( UM const& U, ZM Z ) const {
        assert( U.num_rows()==size() ) ;
        assert( Z.num_rows()==size()-1 ) ;

       // fill( Z(0,glas2::all()), 0.0 ) ;
        for (int i=1; i<size(); ++i)
          Z(i-1,glas2::all()) = weights_(i-1) * U( 0, glas2::all() ) + poles_(i-1) * U( i, glas2::all() ) ;
      } // multiply_M()

    public:
      template <typename ValueType>
      using handle_type = partial_fractions_handle< ValueType, poles_type const&, weights_type const&, size_type > ;

      template <typename ValueType>
      handle_type<ValueType> handle() const { return handle_type<ValueType>( poles_, weights_ ) ; }

    private:
      poles_type const&   poles_ ;
      weights_type const& weights_ ;
  } ; // partial_fractions

} } // namespace CORK::basis4cork

#endif
