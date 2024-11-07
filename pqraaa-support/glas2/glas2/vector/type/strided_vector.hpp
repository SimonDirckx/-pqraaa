//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_strided_vector_hpp
#define glas2_vector_type_strided_vector_hpp

#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/algorithm/ops_assign.hpp>
#include <glas2/vector/type/indirect_vector.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/range_from_end.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/vector/concept/strided_dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <algorithm>
#include <type_traits>

namespace glas2 {

  template <typename T, typename S>
  class strided_vector {
    public:
      typedef T value_type ;
      typedef S size_type ;

    public:
        strided_vector( T* ptr, S n, S stride )
        : ptr_(ptr)
        , size_( n )
        , stride_( stride )
        {}

        // Copy reference !!
        strided_vector( strided_vector const& that )
        : ptr_(that.ptr_)
        , size_( that.size_ )
        , stride_( that.stride_ )
        {}


    public:
      size_type size() const { return size_ ; }
      size_type stride() const { return stride_ ; }
      value_type * const ptr() const {return ptr_ ; }

      value_type const& operator() ( size_type i ) const {
        assert( i>=0 && i<size_ ) ;
        return ptr_[stride_*i] ;
      }

      value_type& operator() ( size_type i ) {
        assert( i>=0 && i<size_ ) ;
        return ptr_[stride_*i] ;
      }

      strided_vector& operator=( strided_vector const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

    public:
      template <typename E>
      strided_vector operator=( E const& that ) {
        return assign( *this, that ) ;
      }

/*      template <typename E>
      strided_vector& operator+=( E const& that ) {
        plus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      strided_vector& operator-=( E const& that ) {
        minus_assign( current_backend(), *this, that ) ;
        return *this ;
      }*/

      template <typename E>
      strided_vector& operator*=( E const& that ) {
        multiplies_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      strided_vector& operator/=( E const& that ) {
        divides_assign( current_backend(), *this, that ) ;
        return *this ;
      }

    public:
      template <typename I>
      typename vector_selection< strided_vector, I >::result_type operator()( I const& s ) {
        return vector_selection< strided_vector, I >::apply( *this, s ) ;
      }

    private:
      T* ptr_ ;
      S  size_ ;
      S  stride_ ;
  } ;


  template <typename T, typename S>
  struct glas_concept< strided_vector<T,S> >
  : StridedDenseVector
  {};

  template <typename T, typename S>
  struct vector_selection< strided_vector<T,S>, all > {
    typedef strided_vector<T,S> result_type ;

     static result_type apply( strided_vector<T,S> v, all ) {
      return v ;
    }
  } ;

  template <typename T, typename S>
  struct vector_selection< strided_vector<T,S>, range > {
    typedef strided_vector<T,slice::size_type> result_type ;
   
    static result_type apply( strided_vector<T,S> v, range r ) {
      assert( r(0) + r.size() <= v.size() ) ;
      return strided_vector<T,range::size_type>( &v(r(0)), r.size(), v.stride() ) ;
    }
  } ;

  template <typename T, typename S>
  struct vector_selection< strided_vector<T,S>, range_from_end > {
    typedef strided_vector<T,slice::size_type> result_type ;
   
    static result_type apply( strided_vector<T,S> v, range_from_end r ) {
      assert( r.begin() <= v.size()-r.from_end() ) ;
      return strided_vector<T,range::size_type>( &v(r.begin()), v.size()-r.begin()-r.from_end(), v.stride() ) ;
    }
  } ;

  template <typename T, typename S>
  struct vector_selection< strided_vector<T,S>, slice > {
    typedef strided_vector<T,slice::size_type> result_type ;
   
    static result_type apply( strided_vector<T,S> v, slice r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return strided_vector<T,slice::size_type>( &v(r(0)), r.size(), v.stride()*r.step() ) ;
    }
  } ;

  template <typename T, typename S, typename S2>
  struct vector_selection< strided_vector<T,S>, S2, typename std::enable_if< is_no_range_or_slice<S2>::value >::type > {
    typedef indirect_vector<strided_vector<T,S>, S2> result_type ;
   
    static result_type apply( strided_vector<T,S> v, S2 r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( v, r ) ;
    }
  } ;

} // namespace glas2


#endif
