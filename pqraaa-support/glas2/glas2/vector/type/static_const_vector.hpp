//  (C) Copyright Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_static_const_vector_hpp
#define glas2_vector_type_static_const_vector_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/scalar/container/static_const_scalar.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/range_from_end.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <type_traits>

namespace glas2 {

  template <typename I, typename T, int N>
  class static_const_vector {
    public:
      static_const_vector( I const& n )
      : size_( n )
      {}

    public:
      typedef I                          size_type ;
      typedef static_const_scalar< T, N> value_type ;

      size_type size() const {
        return size_ ;
      }
      value_type operator() ( size_type i ) const {
        assert( i>=0 && i<size_ ) ;
        return value_type() ;
      }

    private:
      size_type  size_ ;
  } ;

  template <typename I, typename T, int N>
  struct glas_concept< static_const_vector<I,T,N> >
  : DenseVector
  {} ;

  template <typename I, typename T, int N>
  struct vector_selection< static_const_vector<I,T,N>, all > {
    typedef static_const_vector<I,T,N> result_type ;

     static result_type apply( static_const_vector<I,T,N> const& v, all ) {
      return v ;
    }
  } ;

  template <typename I, typename T, int N>
  struct vector_selection< static_const_vector<I,T,N>, range > {
    typedef static_const_vector<I,T,N> result_type ;

    static result_type apply( static_const_vector<I,T,N> v, range r ) {
      assert( r.begin() + r.size() <= v.size() ) ;
      return result_type( r.size() ) ;
    }
  } ;


  template <typename I, typename T, int N>
  struct vector_selection< static_const_vector<I,T,N>, range_from_end > {
    typedef static_const_vector<I,T,N> result_type ;

     static result_type apply( static_const_vector<I,T,N> v, range_from_end r ) {
      assert( r.begin() <= v.size()-r.from_end() ) ;
      return result_type( v.size()-r.begin()-r.from_end() ) ;
    }
  } ;


  template <typename I, typename T, int N>
  struct vector_selection< static_const_vector<I,T,N>, slice > {
    typedef static_const_vector<I,T,N> result_type ;
   
     static result_type apply( static_const_vector<I,T,N> v, slice r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( r.size() ) ;
    }
  } ;

  template <typename I, typename T, int N, typename S2>
  struct vector_selection< static_const_vector<I,T,N>, S2, typename std::enable_if< is_no_range_or_slice<S2>::value >::type > {
    typedef static_const_vector<I,T,N> result_type ;
   
    static result_type apply( static_const_vector<I,T,N> v, S2 r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( r.size() ) ;
    }
  } ;

} // namespace glas2

#endif
