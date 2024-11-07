//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_container_static_matrix_hpp
#define glas2_matrix_container_static_matrix_hpp

#include <glas2/matrix/concept/matrix_transform.hpp>
#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <type_traits>
#include <cassert>
#include <memory>
#ifndef NDEBUG
#include <algorithm>
#include <limits>
#endif

namespace glas2 {

  template <typename T, int M, int N, typename O=column_major>
  class static_matrix
  : public contiguous_matrix< T, int, O >
  {
    public:
      typedef contiguous_matrix< T, int, O > base_type ;
      typedef T              value_type ;
      typedef std::ptrdiff_t size_type ;

    public:
      static_matrix()
      : base_type( ptr_, M, N )
      {
#if !defined(NDEBUG) || defined(GLAS_INITIALIZE_NAN)
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+M*N, std::numeric_limits<T>::quiet_NaN() ) ;
#elif defined(GLAS_INITIALIZE_ZERO)
        std::fill( this->ptr(), this->ptr()+M*N, 0.0 ) ;
#endif
      }

      // Deep copy !!
      static_matrix( static_matrix const& that )
      : base_type( ptr_, M, N )
      {
        std::copy( that.ptr_, that.ptr_+M*N, ptr_ ) ;
      }

    public:
      base_type operator=( static_matrix const& that ) {
        return base_type(*this) = that ;
      }

      template <typename E>
      base_type operator=( E const& that ) {
        return base_type(*this) = that ;
      }

    private:
      value_type ptr_[M*N] ;
  } ;


  template <typename T, int M, int N, typename O>
  struct glas_concept< static_matrix<T,M,N,O> >
  : glas_concept< typename static_matrix<T,M,N,O>::base_type >
  {};

  template <typename Tag, typename T, int M, int N, typename O>
  struct matrix_transform< Tag, static_matrix<T,M,N,O> >
  : matrix_transform< Tag, typename static_matrix<T,M,N,O>::base_type >
  {} ;

  template <typename T, int M, int N, typename O, typename S1, typename S2>
  struct matrix_selection< static_matrix<T,M,N,O>, S1, S2 >
  : matrix_selection< typename static_matrix<T,M,N,O>::base_type, S1, S2 >
  {} ;

} // namespace glas2

#endif
