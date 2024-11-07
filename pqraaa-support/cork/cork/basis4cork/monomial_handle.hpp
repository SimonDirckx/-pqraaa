//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_monomial_handle_hpp
#define cork_basis4cork_monomial_handle_hpp

#include <cork/basis/monomial.hpp>
#include <glas2/matrix.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename T, typename I>
  class monomial_handle
  {
    public:
      typedef T value_type ;
      typedef I size_type ;

    public:
      explicit monomial_handle( size_type size )
      : size_( size )
      {}

    public:
      void shift( value_type s ) { s_ = s ; }

      value_type shift() const { return s_ ; }

      size_type size() const { return size_ ; }

    public:
      value_type phi_0() const { return 1. ; }

    public:
      template <typename Z0, typename ZZ>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z ) const {
        assert( Z.num_rows()==size_-1 ) ;
        fill( Z, 0.0 ) ;
        if (Z.num_rows()>0) Z( 0, glas2::all() ) = -s_ * z0 ;
      } // lower_solve_right_hand_side()

      template <typename ZM, typename Backend=glas2::default_backend>
      void solve( ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assert( Z.num_rows()==size_-1 ) ;
        for (typename ZM::size_type i=1; i<Z.num_rows(); ++i) {
          auto Z1 = Z(i, glas2::all()) ;
          plus_assign( backend, Z1, s_ * Z(i-1, glas2::all()) ) ;
        }
      } // solve()

    private:
      size_type  size_ ;
      value_type s_ ;
  } ; // monomial

} } // namespace CORK::basis4CORK

#endif
