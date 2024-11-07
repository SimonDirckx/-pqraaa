//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_container_dense_vector_hpp
#define glas2_vector_container_dense_vector_hpp

#include <glas2/vector/type/contiguous_vector.hpp>
#include <type_traits>
#include <cassert>
#ifndef NDEBUG
#include <algorithm>
#include <limits>
#endif

namespace glas2 {

  template <typename T>
  class dense_vector
  : public contiguous_vector< T, std::ptrdiff_t >
  {
    public:
      typedef T                                      value_type ;
      typedef std::ptrdiff_t                         size_type ;
      typedef contiguous_vector< T, std::ptrdiff_t > base_type ;

    public:
      dense_vector( size_type n )
      : v_( new T[n] )
      , owner_(true)
      {
        this->reset( v_, n ) ;
#ifndef NDEBUG
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( v_, v_+n, std::numeric_limits<T>::quiet_NaN() ) ;
#endif
      }

      ~dense_vector() {
        if (owner_) delete [] v_ ;
      }

    private:
      // Do not allow copy constructor
      dense_vector( dense_vector const& that ) = delete ;

    public:
      // Copies of the vector become invalid
      dense_vector& resize( size_type n ) {
        assert( owner_ ) ;
        delete [] v_ ;
        v_ = new T[n] ;
#ifndef NDEBUG
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( v_, v_+n, std::numeric_limits<T>::quiet_NaN() ) ;
#endif
        this->reset( v_, n ) ;
        return *this ;
      }

    public:
      contiguous_vector< T, std::ptrdiff_t > reference() { return base_type(*this) ; }

    public:
      base_type& operator=( dense_vector const& that ) {
        base_type(*this) = that ;
        return *this ;
      }

      base_type& operator=( dense_vector&& that ) = delete ;

    public:
      template <typename E>
      base_type& operator=( E const& that ) {
        return base_type(*this) = that ;
      }

    private:
      T*   v_ ;
  } ;


  template <typename T>
  struct glas_concept< dense_vector<T> >
  : glas_concept< typename dense_vector<T>::base_type >
  {};

} // namespace glas2

#endif
