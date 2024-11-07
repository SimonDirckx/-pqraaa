//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_dense_hpp
#define cork_coefficient_matrices_dense_hpp

#include <cork/utility/ref.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cork/concept/has_type_member.hpp>
#include <cork/backend/default_backend.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  // Square dense matrices.
  template <typename Matrices, typename Backend=glas2::default_backend>
  class dense
  {
    public:
      typedef typename CORK::deref_type< Matrices>::type             matrices_type ;
      typedef typename matrices_type::value_type      value_type ;
      typedef typename matrices_type::size_type       grade_type ;
      typedef typename matrices_type::size_type       size_type ;

      template <typename T>
      using value_type_for = typename std::common_type< value_type, T >::type ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< value_type, T > > ;

    public:
      dense( Matrices A, Backend const& backend=glas2::default_backend() )
      : A_( A )
      , backend_( backend )
      {
        assert( A_.num_rows() != 0 ) ;
        assert( A_.num_columns() % A_.num_rows() == 0 ) ;
      }

    public:
      size_type num_rows() const { return A_.num_rows() ; }
      size_type num_columns() const { return A_.num_rows() ; }

      grade_type num_matrices() const { return A_.num_columns()/A_.num_rows() ; }

    public:
      matrices_type const& matrices() const { return A_ ; }

      auto operator() ( int i ) const {
        return A_( glas2::all(), glas2::range(i*num_columns(), i*num_columns()+num_columns()) ) ;
      }

      template <typename M>
      void fill( int i, M m ) const {
        m = A_( glas2::all(), glas2::range(i*num_columns(), i*num_columns()+num_columns()) ) ;
      }

    public:
      template <typename X, typename W>
      void multiply_add( size_type i, X const& x, W w ) const {
        assert( i>=0 && i< num_matrices() ) ;
        //w += glas2::multiply( A_( glas2::all(), glas2::range( i*num_columns(), (i+1)*num_columns() ) ), x ) ;
        plus_assign( backend_, w, glas2::multiply( A_( glas2::all(), glas2::range( i*num_columns(), (i+1)*num_columns() ) ), x ) ) ;
      } // apply_scheduled()

      template <typename Coefs, typename Matrix>
      void accumulate( Coefs const& coefs, Matrix A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
        assert( coefs.size()==num_matrices() ) ;

        for ( int i=0; i<coefs.size(); ++i ) {
           //A += coefs(i) * A_( glas2::all(), glas2::range(i*num_rows(),(i+1)*num_rows()) ) ;
           plus_assign( backend_, A, coefs(i) * A_( glas2::all(), glas2::range(i*num_rows(),(i+1)*num_rows()) ) ) ;
        }
      } // accumulate()

    private:
      Matrices       A_ ;
      Backend const& backend_ ;
  } ; // dense

  template <typename Matrices>
  dense<Matrices> const& make_dense( dense<Matrices> const& d ) { return d; }

} } // namespace CORK::coefficient_matrices

#endif
