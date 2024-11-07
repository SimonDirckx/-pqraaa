//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_indirect_vector_hpp
#define glas2_vector_type_indirect_vector_hpp

#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/algorithm/ops_assign.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/algorithm/assign.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <algorithm>
#include <type_traits>

namespace glas2 {

  template <typename V, typename S>
  class indirect_vector {
    static_assert( is<DenseVector,S>::value, "S must be a vector" );

    public:
      typedef V                      vector_type ;
      typedef S                      selection_type ;
      typedef typename V::value_type value_type ;
      typedef typename S::size_type  size_type ;

    public:
        indirect_vector( V vector, S selection )
        : vector_( vector )
        , selection_( selection )
        {}

        // Copy reference !!
        indirect_vector( indirect_vector const& that )
        : vector_( that.vector_ )
        , selection_( that.selection_ )
        {}


    public:
      size_type size() const { return selection_.size() ; }
      size_type selection() const { return selection_ ; }

      value_type const& operator() ( size_type i ) const {
        assert( i>=0 && i<size() ) ;
        return vector_( selection_(i) ) ;
      }

      value_type& operator() ( size_type i ) {
        assert( i>=0 && i<size() ) ;
        return vector_( selection_(i) ) ;
      }

      indirect_vector& operator=( indirect_vector const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

    public:
      template <typename E>
      indirect_vector operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

      template <typename E>
      indirect_vector& operator+=( E const& that ) {
        plus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      indirect_vector& operator-=( E const& that ) {
        minus_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      indirect_vector& operator*=( E const& that ) {
        multiplies_assign( current_backend(), *this, that ) ;
        return *this ;
      }

      template <typename E>
      indirect_vector& operator/=( E const& that ) {
        divides_assign( current_backend(), *this, that ) ;
        return *this ;
      }

    public:
      template <typename I>
      typename vector_selection< indirect_vector, I >::result_type operator()( I const& s ) {
        return vector_selection< indirect_vector, I >::apply( *this, s ) ;
      }

    private:
      vector_type    vector_ ;
      selection_type selection_ ;
  } ;


  template <typename V, typename S>
  struct glas_concept< indirect_vector<V,S> >
  : DenseVector
  {};


  template <typename V, typename S>
  struct vector_selection< indirect_vector<V,S>, range > {
    typedef typename vector_selection< S, range >::result_type selected_range_type ;
    typedef indirect_vector<V, selected_range_type> result_type ;
   
    static result_type apply( indirect_vector<V,S> v, range r ) {
      assert( r(0) + r.size() <= v.size() ) ;
      return result_type( v.vector(), v.selection( r ) ) ;
    }
  } ;

  template <typename V, typename S>
  struct vector_selection< indirect_vector<V,S>, slice > {
    typedef typename vector_selection< S, slice >::result_type selected_slice_type ;
    typedef indirect_vector<V, selected_slice_type> result_type ;
   
    static result_type apply( indirect_vector<V,S> v, slice r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( v.vector(), v.selection( r ) ) ;
    }
  } ;

  template <typename V, typename S, typename S2>
  struct vector_selection< indirect_vector<V,S>, S2, typename std::enable_if< is_no_range_or_slice<S2>::value >::type > {
    typedef typename vector_selection< S, S2 >::result_type selected_indirect_type ;
    typedef indirect_vector<V, selected_indirect_type> result_type ;
   
    static result_type apply( indirect_vector<V,S> v, S2 r ) {
      assert( r(r.size()-1) < v.size() ) ;
      return result_type( v.vector(), v.selection( r ) ) ;
    }
  } ;

} // namespace glas2


#endif
