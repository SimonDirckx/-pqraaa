//  C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2019.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_rational_newton_handle_hpp
#define cork_basis4cork_rational_newton_handle_hpp

#include <cork/basis/rational_newton.hpp>
#include <cork/exception/linear_solver_failure.hpp>
#include <glas2/matrix.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#include <limits>
#include <cmath>

namespace CORK { namespace basis4CORK {

  template <typename T, typename Points, typename I>
  class rational_newton_handle
  {
    public:
      typedef typename std::decay< Points >::type                                   points_type ;
      typedef typename std::common_type< typename points_type::value_type, T>::type value_type ;
      typedef decltype(std::abs(value_type()))                                      real_type ;
      typedef I                                                                     size_type ;

    public:
      explicit rational_newton_handle( Points nodes, Points poles, Points scaling )
      : nodes_( nodes )
      , poles_( poles )
      , scalings_( scaling )
      {}

    public:
      void shift( value_type s ) {
        if (prod(s-poles_)==0.0) throw exception::linear_solver_failure() ;
        s_ = s ;
      }

      value_type shift() const { return s_ ; }

      size_type size() const { return nodes_.size()-1 ; }

    public:
      value_type phi_0() const {
        value_type factor = 1./scalings_(0) ;
        if (poles_(size()-1)!=std::numeric_limits<real_type>::infinity()) factor /= (poles_(size()-1)-s_) ;
        return factor ;
      }

    public:
      template <typename Z0, typename ZZ>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z ) const {
        assert( Z.num_rows()==size()-1 ) ;
        fill( Z, 0.0 ) ;
        Z( 0, glas2::all() ) = (nodes_(0)-s_) * z0 ;
      } // lower_solve_right_hand_side()

      template <typename ZM>
      void solve( ZM Z ) const {
        assert( Z.num_rows()==size()-1 ) ;
        Z(0, glas2::all() ) /= scalings_(1)*(poles_(0)-s_) ;
        for (typename ZM::size_type i=1; i<Z.num_rows(); ++i) {
          Z(i, glas2::all() ) += (s_-nodes_(i)) * Z(i-1, glas2::all() ) ;
          Z(i, glas2::all() ) /= scalings_(i+1)*(poles_(i)-s_) ;
        }
      } // solve()

    private:
      Points     nodes_ ;
      Points     poles_ ;
      Points     scalings_ ;
      value_type s_ ;
  } ; // rational_newton_handle

} } // namespace CORK::basis4CORK

#endif
