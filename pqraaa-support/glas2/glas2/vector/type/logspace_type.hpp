//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_type_logspace_type_hpp
#define glas2_vector_type_logspace_type_hpp

#include <glas2/vector/type/linspace_type.hpp>

namespace glas2 {

  template <typename T, typename S=int>
  class logspace_type
  : public linspace_type<T,S>
  {
    public:
      typedef T                                value_type ;
      typedef S                                size_type ;

    public:
      inline logspace_type( value_type const& begin, value_type const& end, size_type number )
      : linspace_type<T,S>( begin, end, number )
      {}

    public:
      inline value_type operator()( size_type i ) const {
        return std::pow(10., static_cast<linspace_type<T,S>const&>(*this)( i ) ) ;
      }
  } ;

  template <typename T, typename S>
  struct glas_concept< logspace_type<T,S> >
  : glas_concept< linspace_type<T,S> >
  {} ;

} // namespace glas2

#endif
