//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_container_shared_vector_hpp
#define glas2_vector_container_shared_vector_hpp

#include <glas2/algorithm/copy.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
#include <type_traits>
#include <memory>
#include <cassert>
#ifndef NDEBUG
#include <algorithm>
#include <limits>
#endif

namespace glas2 {

  template <typename T, typename S=std::ptrdiff_t>
  class shared_vector
  : public contiguous_vector< T, S >
  {
    public:
      typedef T                         value_type ;
      typedef S                         size_type ;
      typedef contiguous_vector< T, S > base_type ;

    public:
      shared_vector( size_type n )
      : base_type( new T[n], n )
      , data_ptr_( this->ptr(), [](value_type* p) {delete[] p; } )
      {
#if !defined(NDEBUG) || defined(GLAS_INITIALIZE_NAN)
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+n, std::numeric_limits<T>::quiet_NaN() ) ;
#elif defined(GLAS_INITIALIZE_ZERO)
        std::fill( this->ptr(), this->ptr()+n, 0.0 ) ;
#endif
      }

      // Is a deep copy!!
      shared_vector( std::initializer_list<value_type> that )
      : base_type( new T[that.size()], that.size() )
      , data_ptr_( this->ptr(), [](value_type* p) {delete[] p; })
      {
        std::copy( that.begin(), that.end(), this->ptr() ) ;
      }

      // Copy reference !!
      shared_vector( shared_vector const& that )
      : base_type( that )
      , data_ptr_( that.data_ptr_ )
      {}

      shared_vector( shared_vector& that )
      : base_type( that )
      , data_ptr_( that.data_ptr_ )
      {}

      // Deep copy
      template <typename E>
      shared_vector( copy_of<E> const& that )
      : base_type( new T[that.expression().size()], that.expression().size() )
      , data_ptr_( this->ptr(), [](value_type* p) {delete[] p; })
      {
        *this = that.expression() ;
      }

      // Move constructor
/*      shared_vector( shared_vector&& that )
      : base_type( that )
      , owner_(true)
      {
        that.owner_ = false ;
      }
*/
    public:
      // Copies of the vector become invalid
      shared_vector& resize( size_type n ) {
        assert( data_ptr_.unique() ) ;
        data_ptr_.reset( new T[n], [](value_type* p) {delete[] p; } ) ;
        this->reset( &(*data_ptr_), n ) ;
#ifndef NDEBUG
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+n, std::numeric_limits<T>::quiet_NaN() ) ;
#endif
        return *this ;
      }

      base_type pass_ref() const { return *this ; }

    public:
      base_type& operator=( shared_vector const& that ) {
        base_type(*this) = that ;
        return *this ;
      }

      base_type operator=( shared_vector&& that ) {
        data_ptr_ = that.data_ptr_ ;
        this->reset( &(*data_ptr_), that.size() ) ;
        return *this ;
      }

      // Deep copy
      shared_vector copy() const {
        shared_vector that( this->size() ) ;
        that = *this ;
        return that ;
      }

    public:
      template <typename E>
      base_type operator=( E const& that ) {
        return base_type(*this) = that ;
      }

    private:
      std::shared_ptr<value_type> data_ptr_ ;
  } ;


  template <typename T, typename S>
  struct glas_concept< shared_vector<T,S> >
  : glas_concept< typename shared_vector<T,S>::base_type >
  {};

  template <typename T, typename S, typename S2>
  struct vector_selection< shared_vector<T,S>, S2 >
  : vector_selection< typename shared_vector<T,S>::base_type, S2 >
  {} ;

} // namespace glas2

#endif
