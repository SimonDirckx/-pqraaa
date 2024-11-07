//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_newton_hpp
#define cork_basis4cork_newton_hpp

#include <cork/basis/newton.hpp>
#include <cork/basis4cork/newton_handle.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename Points, typename I>
  class basis4CORK< basis::newton<Points,I> >
  {
    public:
      typedef typename std::decay<Points>::type points_type ;
      typedef I                                 size_type ;

      template< typename T>
      using handle_type = newton_handle<T, points_type const&, size_type> ;

      template< typename T>
      using value_type = typename handle_type<T>::value_type ;

    public:
      explicit basis4CORK( basis::newton<Points,I> const& basis )
      : points_( basis.points() )
      {}

    public:
      // Basis4CORK
      int num_terms() const { return points_.size()+1 ; }
      int size() const { return points_.size() ; }

      typename std::add_lvalue_reference< typename std::add_const<Points>::type >::type points() const& {
        return points_ ;
      }

    public:
      template <typename T>
      handle_type<T> handle() const {
        return handle_type<T>( points_ ) ;
      }

    public:
      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_M( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        fill( A, 0.0 ) ;
        fill( diagonal( A, 1 ), 1.0 ) ;
        diagonal( A, 0 ) = points_(glas2::range(0,size()-1)) ;
      } // fill_M

      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, AA >::value >::type fill_M( AA& A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        glas2::push_back( A, speye(size()-1,size()-1), glas2::range(0,size()-1), glas2::range(1,size()) ) ;
        for (int i=0; i<size()-1; ++i) A.push_back( i, i, points_(i) ) ;
      } // fill_M

      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_N( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        fill( A, 0.0 ) ;
        fill( diagonal( A, 0 ), 1.0 ) ;
      } // fill_N

      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, AA >::value >::type fill_N( AA& A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        glas2::push_back( A, speye(size()-1,size()-1), glas2::range(0,size()-1), glas2::range(0,size()-1) ) ;
      } // fill_N

    public:
      // Basis4CORK
      template <typename UM, typename ZM>
      void multiply_N( UM const& U, ZM Z ) const {
        Z = U( glas2::range(0,U.num_rows()-1), glas2::all() ) ;
      } // multiply_N()

      // Basis4CORK
      template <typename UM, typename ZM>
      void multiply_M( UM const& U, ZM Z ) const {
        for (size_type i=0; i<Z.num_rows(); ++i)
          Z(i, glas2::all()) = points_(i) * U(i, glas2::all()) ;
        Z += U(glas2::range_from_end(1,0), glas2::all()) ;
      } // multiply_M()

    private:
      points_type const& points_ ;
  } ; // newton

} } // namespace CORK::basis4cork

#endif
