//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_unit_vector_hpp
#define glas2_vector_type_unit_vector_hpp

#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <algorithm>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T>
  class unit_vector {
    public:
      typedef T   value_type ;
      typedef int size_type ;

    public:
        unit_vector()
        : size_(0)
        , pos_(0)
        {}

        unit_vector( size_type n, size_type pos )
        : size_( n )
        , pos_( pos )
        {}

        // Copy reference !!
        unit_vector( unit_vector const& that )
        : size_( that.size_ )
        , pos_( that.pos_ )
        {}

    public:
        void reset( size_type n, size_type pos ) {
          size_ = n ;
          pos_ = pos ;
        }

    public:
      size_type size() const { return size_ ; }
      size_type position() const { return pos_ ; }

      value_type operator() ( size_type i ) const {
        assert( i>=0 && i<size_ ) ;
        if (i==pos_) return 1 ;
        return 0 ;
      }

      template <typename I>
      typename vector_selection< unit_vector, I >::result_type operator()( I const& s ) {
        return vector_selection< unit_vector, I >::apply( *this, s ) ;
      }

    private:
      size_type  size_ ;
      size_type  pos_ ;
  } ;


  template <typename T>
  struct glas_concept< unit_vector<T> >
  : DenseVector
  {};


  template <typename T>
  struct vector_selection< unit_vector<T>, all > {
    typedef unit_vector<T> result_type ;

     static result_type apply( unit_vector<T> v, all ) {
       return v ;
    }
  } ;


  template <typename T>
  struct vector_selection< unit_vector<T>, range > {
    typedef unit_vector<T> result_type ;

     static result_type apply( unit_vector<T> v, range r ) {
      assert( r(0) + r.size() <= v.size() ) ;
      return result_type( r.size(), v.position()-r.begin() ) ;
    }
  } ;


  template <typename T>
  struct vector_selection< unit_vector<T>, slice > {
    typedef unit_vector<T> result_type ;
   
     static result_type apply( unit_vector<T> v, slice r ) {
      assert( r(r.size()-1) < v.size() ) ;
      int pos = v.position() - r.begin() ;
      if (pos%r.step()==0) pos = pos / r.step() ;
      else pos = -1 ;
      return result_type( r.size(), pos ) ;
    }
  } ;

} // namespace glas2


#endif
