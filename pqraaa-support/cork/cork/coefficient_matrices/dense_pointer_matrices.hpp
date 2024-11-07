//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_dense_pointer_matrices_hpp
#define cork_coefficient_matrices_dense_pointer_matrices_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cork/backend/default_backend.hpp>
#include <cstddef>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  // Square dense matrices.
  template <typename ValueType> 
  class dense_pointer_matrices
  {
    public:
      typedef ValueType      value_type ;
      typedef ptrdiff_t      grade_type ;
      typedef ptrdiff_t      size_type ;

    public:
      dense_pointer_matrices( ValueType *matsptr, grade_type num_mats, size_type mats_offset, size_type n, size_type lda )
      : matsptr_( matsptr )
      , num_mats_( num_mats )
      , mats_offset_( mats_offset )
      , n_( n )
      , lda_( lda )
      {
      }

    public:
      size_type num_rows() const { return n_ ; }
      size_type num_columns() const { return n_ ; }

      grade_type num_matrices() const { return num_mats_ ; }

    public:

//      auto operator() ( int i ) const {
//        return A_( glas2::all(), glas2::range(i*num_columns(), i*num_columns()+num_columns()) ) ;
//      }

//      template <typename M>
//      void fill( int i, M m ) const {
//        m = A_( glas2::all(), glas2::range(i*num_columns(), i*num_columns()+num_columns()) ) ;
//      }

    public:
//      template <typename X, typename W>
//      void multiply_add( size_type i, X const& x, W& w ) const {
//        assert( i>=0 && i< num_matrices() ) ;
//        //w += glas2::multiply( A_( glas2::all(), glas2::range( i*num_columns(), (i+1)*num_columns() ) ), x ) ;
//        plus_assign( backend_, w, glas2::multiply( A_( glas2::all(), glas2::range( i*num_columns(), (i+1)*num_columns() ) ), x ) ) ;
//      } // apply_scheduled()
      template <typename X, typename W>
      void multiply_add( size_type i, X const& x, W& w ) const {
        for (size_type c=0; c<n_; c++) {
          for (size_type r=0; r<n_; r++) { w(r) += x(c)* (matsptr_[i*mats_offset_+c*lda_+r]); }
        }
      } // multiply_add()

      template <typename Shift, typename X>
      void solve( Shift const& shift, X x, bool is_new_shift ) {
//        static_assert( glas2::is<glas2::ContiguousDenseVector,X>::value, "X should be a contiguous vector" ) ;
//        assert( x.size()==num_rows() ) ;
//        solve_( shift, x.size(), boost::numeric::bindings::begin_value(x), is_new_shift ) ;
      } // solve()

//      template <typename Coefs, typename Matrix>
//      void accumulate( Coefs const& coefs, Matrix& A ) const {
//        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
//        assert( coefs.size()==num_matrices() ) ;
//
//        for ( int i=0; i<coefs.size(); ++i ) {
//           //A += coefs(i) * A_( glas2::all(), glas2::range(i*num_rows(),(i+1)*num_rows()) ) ;
//           plus_assign( backend_, A, coefs(i) * A_( glas2::all(), glas2::range(i*num_rows(),(i+1)*num_rows()) ) ) ;
//        }
//      } // accumulate()

    private:
      ValueType* matsptr_;
      grade_type num_mats_;
      size_type  mats_offset_;
      size_type  n_;
      size_type  lda_;
  } ; // dense

  template <typename ValueType>
  dense_pointer_matrices<ValueType> const& make_dense_pointer( dense_pointer_matrices<ValueType> const& d ) { return d; }

} } // namespace CORK::coefficient_matrices

#endif
