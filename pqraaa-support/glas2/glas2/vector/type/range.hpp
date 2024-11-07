//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_range_hpp
#define glas2_vector_type_range_hpp

#include <glas2/vector/type/all.hpp>
#include <glas2/vector/algorithm/vector_selection.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/concept/is_no_range_or_slice.hpp>
#include <glas2/concept/concept.hpp>

namespace glas2 {

  class range {
    public:
      typedef int value_type ;
      typedef int size_type ;

    public:
      inline range( value_type begin, value_type end )
      : begin_( begin )
      , end_( end )
      {}

    public:
      inline size_type begin() const { return begin_ ; }
      inline size_type end() const { return end_ ; }
      inline size_type size() const { return end_-begin_ ; }

      inline value_type operator()( size_type i ) const { return begin_ + i ; }

      template <typename I>
      typename std::enable_if< !std::is_integral<I>::value, typename vector_selection< range, I >::result_type>::type operator()( I const& s ) const {
        return vector_selection< range, I >::apply( *this, s ) ;
      }

    private:
      value_type begin_ ;
      value_type end_ ;
  } ;

  template <>
  struct glas_concept< range >
  : DenseVector
  {} ;

  template <>
  struct is_no_range_or_slice< range >
  : std::false_type
  {} ;

  template <>
  struct vector_selection< range, all > {
    typedef range result_type ;
   
    static result_type apply( range v, all ) {
      return v ;
    }
  } ;

  template <>
  struct vector_selection< range, range > {
    typedef range result_type ;
   
    static result_type apply( range v, range r ) {
      assert( (r.size()==0) || (r(r.size()-1) < v.size()) ) ;
      return range( v.begin()+r.begin(), v.begin()+r.end() ) ;
    }
  } ;


} // namespace glas2

#endif
