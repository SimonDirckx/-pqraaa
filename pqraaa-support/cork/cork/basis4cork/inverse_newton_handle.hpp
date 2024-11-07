//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_newton_handle_hpp
#define cork_basis4cork_newton_handle_hpp

#include <cork/basis/newton.hpp>
#include <glas2/matrix.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename T, typename Points, typename I>
  class newton_handle
  {
    public:
      typedef typename std::decay< Points >::type                                   points_type ;
      typedef typename std::common_type< typename points_type::value_type, T>::type value_type ;
      typedef I                                                                     size_type ;

    public:
      explicit newton_handle( Points points )
      : points_( points )
      {}

    public:
      void shift( value_type s ) { s_ = s ; }

      value_type shift() const { return s_ ; }

      size_type size() const { return points_.size() ; }

    public:
      template <typename Z0, typename ZZ>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z ) const {
        assert( Z.num_rows()==size()-1 ) ;
        fill( Z, 0.0 ) ;
        Z( 0, glas2::all() ) = (points_(0)-s_) * z0 ;
      } // lower_solve_right_hand_side()

      template <typename ZM>
      void solve( ZM Z ) const {
        assert( Z.num_rows()==size()-1 ) ;
        for (typename ZM::size_type i=1; i<Z.num_rows(); ++i) {
          Z(i, glas2::all() ) += (s_-points_(i)) * Z(i-1, glas2::all() ) ;
        }
      } // solve()

    private:
      Points     points_ ;
      value_type s_ ;
  } ; // newton_handle

} } // namespace CORK::basis4CORK

#endif
