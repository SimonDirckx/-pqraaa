//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2019.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_rational_newton_hpp
#define cork_basis4cork_rational_newton_hpp

#include <cork/basis/rational_newton.hpp>
#include <cork/basis4cork/rational_newton_handle.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename Points>
  class basis4CORK< basis::rational_newton<Points> >
  {
    public:
      typedef typename std::decay<Points>::type points_type ;
      typedef typename points_type::size_type   size_type ;

      template< typename T>
      using handle_type = rational_newton_handle<T, points_type const&, size_type> ;

      template< typename T>
      using value_type_for = typename handle_type<T>::value_type ;

    public:
      explicit basis4CORK( basis::rational_newton<Points> const& basis )
      : nodes_( basis.nodes() )
      , poles_( basis.poles() )
      , scaling_( basis.scaling() )
      {}

    public:
      // Basis4CORK
      int num_terms() const { return nodes_.size() ; }
      int size() const { return nodes_.size()-1 ; }

      typename std::add_lvalue_reference< typename std::add_const<Points>::type >::type nodes() const& {
        return nodes_ ;
      }

      typename std::add_lvalue_reference< typename std::add_const<Points>::type >::type poles() const& {
        return poles_ ;
      }

      typename std::add_lvalue_reference< typename std::add_const<Points>::type >::type scaling() const& {
        return scaling_ ;
      }

    public:
      template <typename T>
      handle_type<T> handle() const {
        return handle_type<T>( nodes_, poles_, scaling_ ) ;
      }

    public:
      // Basis4CORK
      template <typename UM, typename ZM>
      void multiply_N( UM const& U, ZM Z ) const {
        assert( Z.num_columns() == U.num_columns() ) ;
        assert( Z.num_rows() == U.num_rows()-1 ) ;
        Z = U( glas2::range_from_end(0,1), glas2::all() ) ;
        for (size_type i=0; i<Z.num_rows(); ++i)
          Z(i,glas2::all()) += scaling_(i+1) * U( i+1, glas2::all() ) ;
      } // multiply_N()

      // Basis4CORK
      template <typename UM, typename ZM>
      void multiply_M( UM const& U, ZM Z ) const {
        assert( Z.num_columns() == U.num_columns() ) ;
        assert( Z.num_rows() == U.num_rows()-1 ) ;
        for (size_type i=0; i<Z.num_rows(); ++i)
          Z(i, glas2::all()) = scaling_(i+1)*poles_(i) * U(i+1, glas2::all()) + nodes_(i) * U(i,glas2::all()) ;
      } // multiply_M()

    private:
      points_type const& nodes_ ;
      points_type const& poles_ ;
      points_type const& scaling_ ;
  } ; // rational_newton

} } // namespace CORK::basis4cork

#endif
