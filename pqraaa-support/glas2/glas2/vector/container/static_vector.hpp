//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_container_static_vector_hpp
#define glas2_vector_container_static_vector_hpp

#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/vector/algorithm/assign.hpp>
#include <type_traits>
#include <cassert>
#ifndef NDEBUG
#include <algorithm>
#include <limits>
#endif

namespace glas2 {

  template <typename T, int N>
  class static_vector
  : public contiguous_vector< T, int >
  {
    public:
      typedef T                           value_type ;
      typedef int                         size_type ;
      typedef contiguous_vector< T, int > base_type ;

    public:
      static_vector( int n=N )
      {
        assert( n==N ) ;
        this->reset( v_, N ) ;
#if !defined(NDEBUG) || defined(GLAS_INITIALIZE_NAN)
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+n, std::numeric_limits<T>::quiet_NaN() ) ;
#elif defined(GLAS_INITIALIZE_ZERO)
        std::fill( this->ptr(), this->ptr()+n, 0.0 ) ;
#endif
      }

      // Is a deep copy!!
      static_vector( std::initializer_list<value_type> that )
      {
        assert( that.size() == N ) ;
        this->reset( v_, N ) ;
        std::copy( that.begin(), that.end(), v_ ) ;
      }

      // Is a deep copy!!
      static_vector( static_vector const& that )
      {
        std::copy( that.v_, that.v_+N, v_ ) ;
        this->reset( v_, N ) ;
      }

      base_type pass_ref() const { return *this ; }

    public:
      static_vector& operator=( static_vector const& that ) {
        std::copy( that.v_, that.v_+N, v_ ) ;
        return *this ;
      }

    public:
      template <typename E>
      static_vector& operator=( E const& that ) {
        static_cast<base_type&>(*this) = that ;
        //assign( static_cast<base_type&>(*this), that ) ;
        return *this ;
      }

    private:
      T   v_[N] ;
  } ;


  template <typename T, int N>
  struct glas_concept< static_vector<T,N> >
  : ContiguousDenseVector
  {};

} // namespace glas2

#endif
