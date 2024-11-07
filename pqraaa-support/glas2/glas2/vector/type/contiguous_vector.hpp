//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_contiguous_vector_hpp
#define glas2_vector_type_contiguous_vector_hpp

#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/algorithm/ops_assign.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/range_from_end.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/indirect_vector.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <algorithm>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T, typename S>
  class contiguous_vector {
    public:
      typedef T value_type ;
      typedef S size_type ;

    public:
        contiguous_vector()
        : ptr_(0)
        , size_(0)
        {}

        contiguous_vector( T* ptr, S n )
        : ptr_(ptr)
        , size_( n )
        {}

        // Copy reference !!
        contiguous_vector( contiguous_vector const& that )
        : ptr_(that.ptr_)
        , size_( that.size_ )
        {}

    public:
        void reset( T* __restrict__ ptr, S n ) {
          ptr_ = ptr ;
          size_ = n ;
        }

    public:
      size_type size() const { return size_ ; }
      value_type* const __restrict__ ptr() const { return ptr_ ; }

    public: // STL
      value_type const* const __restrict__ begin() const { return ptr_ ; }
      value_type const* const __restrict__ end() const { return ptr_+size_ ; }

      value_type* const __restrict__ begin() { return ptr_ ; }
      value_type* const __restrict__ end() { return ptr_+size_ ; }

    public:
      template <typename I>
      typename std::enable_if< std::is_integral<I>::value, value_type& >::type operator() ( I i ) const {
        assert( size_type(i)>=0 && size_type(i)<size() ) ;
        return ptr_[size_type(i)] ;
      }

      template <typename I>
      typename std::enable_if< std::is_integral<I>::value, value_type& >::type operator[] ( I i ) const {
        assert( size_type(i)>=0 && size_type(i)<size() ) ;
        return ptr_[size_type(i)] ;
      }

      template <typename I>
      typename std::enable_if< std::is_integral<I>::value, value_type& >::type operator() ( I i ) {
        assert( i>=0 && i<size() ) ;
        return ptr_[i] ;
      }

      template <typename I>
      typename std::enable_if< std::is_integral<I>::value, value_type& >::type operator[] ( I i ) {
        assert( i>=0 && i<size() ) ;
        return ptr_[i] ;
      }

      contiguous_vector& operator=( contiguous_vector const& that ) {
        assert( that.size()==size() ) ;
        std::copy( that.ptr_, that.ptr_+that.size(), ptr_ ) ;
        return *this ;
      }

    public:
      template <typename E>
      contiguous_vector& operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

/*      template <typename E>
      contiguous_vector& operator+=( E const& that ) {
        plus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      contiguous_vector& operator-=( E const& that ) {
        minus_assign( current_backend(), *this, that ) ;
        return *this ;
      }*/

      template <typename E>
      contiguous_vector& operator*=( E const& that ) {
        multiplies_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      contiguous_vector& operator/=( E const& that ) {
        divides_assign( current_backend(), *this, that ) ;
        return *this ;
      }

    public:
      template <typename I>
      typename std::enable_if< !std::is_integral<I>::value, typename vector_selection< contiguous_vector, I >::result_type>::type operator()( I const& s ) {
        return vector_selection< contiguous_vector, I >::apply( *this, s ) ;
      }

      template <typename I>
      typename std::enable_if< !std::is_integral<I>::value, typename vector_selection< contiguous_vector, I >::result_type>::type operator()( I const& s ) const {
        return vector_selection< contiguous_vector, I >::apply( *this, s ) ;
      }

    private:
      value_type* __restrict__ ptr_ ;
      size_type   size_ ;
  } ;


  template <typename T, typename S>
  struct glas_concept< contiguous_vector<T,S> >
  : ContiguousDenseVector
  {};


  template <typename T, typename S>
  struct vector_selection< contiguous_vector<T,S>, all > {
    typedef contiguous_vector<T,S> result_type ;

     static result_type apply( contiguous_vector<T,S> v, all ) {
      return v ;
    }
  } ;

  template <typename T, typename S>
  struct vector_selection< contiguous_vector<T,S>, range > {
    typedef contiguous_vector<T,S> result_type ;

     static result_type apply( contiguous_vector<T,S> v, range r ) {
      assert( r.begin() + r.size() <= v.size() ) ;
      return result_type( v.ptr()+r.begin(), r.size() ) ;
    }
  } ;


  template <typename T, typename S>
  struct vector_selection< contiguous_vector<T,S>, range_from_end > {
    typedef contiguous_vector<T,S> result_type ;

     static result_type apply( contiguous_vector<T,S> v, range_from_end r ) {
      assert( r.begin() <= v.size()-r.from_end() ) ;
      return result_type( v.ptr()+r.begin(), v.size()-r.begin()-r.from_end() ) ;
    }
  } ;


  template <typename T, typename S>
  struct vector_selection< contiguous_vector<T,S>, slice > {
    typedef strided_vector<T,S> result_type ;
   
     static result_type apply( contiguous_vector<T,S> v, slice r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( v.ptr()+r.begin(), r.size(), r.step() ) ;
    }
  } ;

  template <typename T, typename S, typename S2>
  struct vector_selection< contiguous_vector<T,S>, S2, typename std::enable_if< is_no_range_or_slice<S2>::value >::type > {
    typedef indirect_vector<contiguous_vector<T,S>, S2> result_type ;
   
    static result_type apply( contiguous_vector<T,S> v, S2 r ) {
      assert( (r.size()==0) || (r(r.size()-1) < v.size()) ) ;
      return result_type( v, r ) ;
    }
  } ;

} // namespace glas2

#endif
