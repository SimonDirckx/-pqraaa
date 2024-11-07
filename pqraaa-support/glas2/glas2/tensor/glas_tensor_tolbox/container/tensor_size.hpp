//  (C) Copyright Karl Meerbergen 2012.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_container_tensor_size_hpp
#define glas_toolbox_tensor_container_tensor_size_hpp

#include <glas/config.hpp>
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER
#include <glas/concept/check_include_level.hpp>

#include <glas/toolbox/tensor/concept/round_bracketed.hpp>
#include <glas/concept/square_bracketed.hpp>
#include <boost/static_assert.hpp>
#include <cassert>
#ifdef GLAS_WARN_CONTAINER_COPY
#ifndef NDEBUG
#include <iostream>
#include <string>
#include <limits>
#endif
#endif

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_CONTAINER


namespace glas {

  namespace detail {

    template <int I, int O>
    struct size_at_functor {
      BOOST_STATIC_ASSERT( (I<=O) ) ;

      template <typename S>
      int& operator() ( S& s ) const {
        return size_at_functor<I,O-1>( s.size_rest_ ) ;
      }

      template <typename S>
      int operator() ( S const& s ) const {
        return size_at_functor<I,O-1>( s.size_rest_ ) ;
      }
    } ;

    template <int O>
    struct size_at_functor {
      template <typename S>
      int& operator() ( S& s ) const {
        return s.size_ ;
      }

      template <typename S>
      int operator() ( S const& s ) const {
        return s.size_ ;
      }
    } ;

  } // namespace detail


  template <int Order>
  struct tensor_size
  {
    template <int I>
    int& at() { return detail::size_at_functor<I,Order>() ( *this ) ; }

    template <int I>
    int at() const { return detail::size_at_functor<I,Order>() ( *this ) ; }

    tensor_size<Order-1> size_rest_ ;
    int                  size_ ;
  } ; // tensor_size

  template <1>
  struct tensor_size
  {
    int size_ ;
  } ; // tensor_size


} // Namespace glas
  

#endif
