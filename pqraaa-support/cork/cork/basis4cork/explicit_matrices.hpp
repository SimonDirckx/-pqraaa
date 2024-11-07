//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_explicit_matrices_hpp
#define cork_basis4cork_explicit_matrices_hpp

#include <cork/basis4cork/explicit_matrices_handle.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getrs.hpp>
//#include <boost/numeric/bindings/lapack/computational/gecon.hpp>
#include <cassert>
#include <string>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename T>
  class explicit_matrices
  {
    private:
      typedef glas2::shared_matrix< T > matrix_type ;
    
    public:
      template <typename ValueType>
      using value_type_for = typename std::common_type< ValueType, T >::type ;

      typedef typename matrix_type::size_type size_type ;

    public:
      explicit explicit_matrices( int size )
      : M_( size-1, size )
      , N_( size-1, size )
      {}

    protected:
      auto const& M() const { return M_ ; }
      auto const& N() const { return N_ ; }

      typename matrix_type::base_type M() { return M_ ; }
      typename matrix_type::base_type N() { return N_ ; }

    public:
      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_M( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        A = M_ ;
      } // fill_M

      //Basis4CORK
/*      template <typename AA>
      typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, AA >::value >::type fill_M( AA A ) {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        glas2::push_back( A, glas2::range(0,size()-1), glas2::range(1,size()), speye(size()-1,size()-1) ) ;
      } // fill_M*/

      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_N( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        A = N_ ;
      } // fill_N

/*      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, AA >::value >::type fill_N( AA A ) {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        glas2::push_back( A, glas2::range(0,size()-1), glas2::range(0,size()-1), speye(size()-1,size()-1) ) ;
      } // fill_N*/

    public:
      size_type size() const { return M_.num_columns() ; }

      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_M( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assign( backend, Z, multiply( M_, U ) ) ;
      } // multiply_M()

      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_N( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assign( backend, Z, multiply( N_, U ) ) ;
      } // multiply_N()

    public:
      template <typename ValueType>
      using handle_type = explicit_matrices_handle< ValueType, T > ;

      template <typename ValueType>
      handle_type<ValueType> handle() const {
        return handle_type<ValueType>( M_, N_ ) ;
      }

    private:
      matrix_type      M_ ;
      matrix_type      N_ ;
  } ; // explicit_matrices

} } // namespace CORK::basis4cork

#endif
