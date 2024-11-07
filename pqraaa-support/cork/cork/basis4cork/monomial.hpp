//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_monomial_hpp
#define cork_basis4cork_monomial_hpp

#include <cork/basis/monomial.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/basis4cork/monomial_handle.hpp>
#include <glas2/matrix.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename I>
  class basis4CORK< basis::monomial<I> >
  {
    public:
      typedef I                     size_type ;

      template< typename T>
      using handle_type = monomial_handle<T, size_type> ;

      template< typename T>
      using value_type = T ;

      template< typename T>
      using value_type_for = T ;

    public:
      explicit basis4CORK( basis::monomial<I> const& basis )
      : size_( basis.num_terms()-1 )
      {}

    public:
      // Basis4CORK
      int num_terms() const { return size_+1 ; }
      int size() const { return size_ ; }

    public:
      template <typename T>
      handle_type<T> handle() const {
        return handle_type<T>( size_ ) ;
      }

    public:
      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_M( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        fill( A, 0.0 ) ;
        fill( diagonal( A, 1 ), 1.0 ) ;
      } // fill_M

      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, AA >::value >::type fill_M( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        glas2::push_back( A, glas2::range(0,size()-1), glas2::range(1,size()), speye(size()-1,size()-1) ) ;
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
      typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, AA >::value >::type fill_N( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        glas2::push_back( A, glas2::range(0,size()-1), glas2::range(0,size()-1), speye(size()-1,size()-1) ) ;
      } // fill_N

    public:
      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_N( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assign( backend, Z, U( glas2::range(0,U.num_rows()-1), glas2::all() ) ) ;
      } // multiply_N()

      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_M( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assign( Z, U( glas2::range_from_end(1,0), glas2::all() ) ) ;
      } // multiply_M()

    private:
      size_type size_ ;
  } ; // monomial

} } // namespace CORK::basis4CORK

#endif
