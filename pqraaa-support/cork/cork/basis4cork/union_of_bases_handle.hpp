//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_union_of_bases_handle_hpp
#define cork_basis4cork_union_of_bases_handle_hpp

#include <cork/basis/union_of_bases.hpp>
#include <cork/basis4cork/union_of_bases_handle.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <type_traits>
#include <cassert>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename B1, typename B2>
  class union_of_bases_handle
  {
    public:
      typedef typename B1::value_type                                                           value_type ;
      typedef typename std::common_type< typename B1::size_type, typename B2::size_type>::type  size_type ;

    public:
      explicit union_of_bases_handle( B1 const& handle_1, B2 const& handle_2 )
      : b1_( handle_1 )
      , b2_( handle_2 )
      {}

    public:
      size_type size() const { return std::max<size_type>(b1_.size(), std::max( b2_.size(), b1_.size()+b2_.size() -1 ) ) ; }

      // Basis4CORK
      void shift( value_type s ) { b1_.shift(s) ; b2_.shift(s) ; }

      // Basis4CORK
      value_type shift() const { return b1_.shift() ; }

    public:
      value_type phi_0() const {
        assert( b2_.phi_0()==b1_.phi_0() );
        return b1_.phi_0() ; 
      }

    public:
      // Basis4CORK
      template <typename Z0, typename ZZ, typename Backend=glas2::default_backend>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z, Backend const& backend=glas2::default_backend() ) const {
        b2_.lower_solve_right_hand_side( z0, Z(glas2::range_from_end(std::max(0,b1_.size()-1),0),glas2::all()) ) ;
        if (b1_.size()>0) b1_.lower_solve_right_hand_side( z0, Z(glas2::range(0,b1_.size()-1), glas2::all()) ) ;
      } // multiply_N()

      // Basis4CORK
      template <typename ZM, typename Backend=glas2::default_backend>
      void solve( ZM Z, Backend const& backend=glas2::default_backend() ) const {
        if (b1_.size()>1) b1_.solve( Z( glas2::range(0,b1_.size()-1), glas2::all() ), backend ) ;
        b2_.solve( Z( glas2::range_from_end(std::max(0,b1_.size()-1),0), glas2::all() ), backend ) ;
      } // solve()

    private:
      B1 b1_ ;
      B2 b2_ ;
  } ; // union_of_bases

} } // namespace CORK::basis4CORK

#endif
