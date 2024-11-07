//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_container_shared_matrix_hpp
#define glas2_matrix_container_shared_matrix_hpp

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

  template <typename T, typename O=column_major>
  class shared_matrix
  : public contiguous_matrix< T, std::ptrdiff_t, O >
  {
    public:
      typedef contiguous_matrix< T, std::ptrdiff_t, O > base_type ;
      typedef T              value_type ;
      typedef std::ptrdiff_t size_type ;

    public:
      shared_matrix( size_type m, size_type n )
      : contiguous_matrix< T, std::ptrdiff_t, O >( new T[m*n], m, n )
      , shared_( this->ptr(), [](value_type* p) {delete[] p;} )
      {
#if !defined(NDEBUG) || defined(GLAS_INITIALIZE_NAN)
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+m*n, std::numeric_limits<T>::quiet_NaN() ) ;
#elif defined(GLAS_INITIALIZE_ZERO)
        std::fill( this->ptr(), this->ptr()+m*n, 0.0 ) ;
#endif
      }

      // Copy reference !!
      shared_matrix( shared_matrix const& that )
      : contiguous_matrix< T, std::ptrdiff_t, O >( that )
      , shared_( that.shared_ )
      {}

    public:
      // Copies of the matrix become invalid
 shared_matrix& resize( size_type n, size_type m ) {
        shared_.reset( new T[n*m],[](value_type* p) {delete[] p; } ) ;
        this->reset( &(*shared_), n, m ) ;
#ifndef NDEBUG
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+n*m, std::numeric_limits<T>::quiet_NaN() ) ;
#endif
        return *this ;
      }

      // Content of the remaining matrix stays
      shared_matrix& safe_resize( size_type m, size_type n ) {
        base_type new_matrix( new T[m*n], m, n ) ;
        new_matrix = (*this)( range(0,m), range(0,n) ) ;
        this->reset( new_matrix.ptr(), m, n ) ;
        delete [] this->ptr() ;
        return *this ;
      }

      base_type const& pass_ref() const{ return *this; }

    public:
      base_type operator=( shared_matrix const& that ) {
        return base_type(*this) = that ;
      }

      // Deep copy
      shared_matrix copy() const {
        shared_matrix that( this->num_rows(), this->num_columns() ) ;
        that = *this ;
        return that ;
      }

      template <typename E>
      base_type operator=( E const& that ) {
        return base_type(*this) = that ;
      }

      /*template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< base_type, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) {
        return matrix_selection< base_type, I1, I2 >::apply( *this, s1, s2 ) ;
      }*/

    private:
      std::shared_ptr<value_type> shared_ ;
  } ;


  template <typename T, typename O>
  struct glas_concept< shared_matrix<T,O> >
  : glas_concept< typename shared_matrix<T,O>::base_type >
  {};

  template <typename Tag, typename T, typename O>
  struct matrix_transform< Tag, shared_matrix<T,O> >
  : matrix_transform< Tag, typename shared_matrix<T,O>::base_type >
  {} ;

  template <typename T, typename O, typename S1, typename S2>
  struct matrix_selection< shared_matrix<T,O>, S1, S2 >
  : matrix_selection< typename shared_matrix<T,O>::base_type, S1, S2 >
  {} ;

} // namespace glas2

#endif
