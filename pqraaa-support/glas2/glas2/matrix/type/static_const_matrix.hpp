//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_static_const_matrix_hpp
#define glas2_matrix_type_static_const_matrix_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/container/static_const_scalar.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/range_from_end.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename T, int N>
  class static_const_matrix {
    public:
      static_const_matrix( I const& m, I const& n )
      : num_rows_( m )
      , num_columns_( n )
      {}

    public:
      typedef I                          size_type ;
      typedef static_const_scalar< T, N> value_type ;

      size_type num_rows() const {
        return num_rows_ ;
      }
      size_type num_columns() const {
        return num_columns_ ;
      }
      value_type operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows_ ) ;
        assert( j>=0 && j<num_columns_ ) ;
        return value_type() ;
      }

    private:
      size_type  num_rows_ ;
      size_type  num_columns_ ;
  } ;

  template <typename I, typename T, int N>
  struct glas_concept< static_const_matrix<I,T,N> >
  : DenseMatrix
  {} ;

/*  template <typename I, typename T, int N>
  struct matrix_selection< static_const_matrix<I,T,N>, all > {
    typedef static_const_matrix<I,T,N> result_type ;

     static result_type apply( static_const_matrix<I,T,N> const& v, all ) {
      return v ;
    }
  } ;

  template <typename I, typename T, int N>
  struct matrix_selection< static_const_matrix<I,T,N>, range > {
    typedef static_const_matrix<I,T,N> result_type ;

    static result_type apply( static_const_matrix<I,T,N> v, range r ) {
      assert( r.begin() + r.size() <= v.size() ) ;
      return result_type( r.size() ) ;
    }
  } ;


  template <typename I, typename T, int N>
  struct matrix_selection< static_const_matrix<I,T,N>, range_from_end > {
    typedef static_const_matrix<I,T,N> result_type ;

     static result_type apply( static_const_matrix<I,T,N> v, range_from_end r ) {
      assert( r.begin() <= v.size()-r.from_end() ) ;
      return result_type( v.size()-r.begin()-r.from_end() ) ;
    }
  } ;


  template <typename I, typename T, int N>
  struct matrix_selection< static_const_matrix<I,T,N>, slice > {
    typedef static_const_matrix<I,T,N> result_type ;
   
     static result_type apply( static_const_matrix<I,T,N> v, slice r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( r.size() ) ;
    }
  } ;

  template <typename I, typename T, int N, typename S2>
  struct matrix_selection< static_const_matrix<I,T,N>, S2, typename std::enable_if< is_no_range_or_slice<S2>::value >::type > {
    typedef static_const_matrix<I,T,N> result_type ;
   
    static result_type apply( static_const_matrix<I,T,N> v, S2 r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( r.size() ) ;
    }
  } ;*/

} // namespace glas2

#endif
