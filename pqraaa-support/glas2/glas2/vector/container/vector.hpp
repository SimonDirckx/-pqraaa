//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_vector_container_vector_hpp
#define glas2_vector_container_vector_hpp

#include <glas2/type/pass_reference.hpp>
#include <glas2/algorithm/copy.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
#include <type_traits>
#include <cassert>
#include <memory>
#include <initializer_list>
#ifndef NDEBUG
#include <algorithm>
#include <limits>
#endif
#ifdef GLAS2_MEMORY_WARNINGS
#include <iostream>
#endif

namespace glas2 {

  template <typename T>
  class vector
  : public contiguous_vector< T, std::ptrdiff_t >
  {
    public:
      typedef T                                      value_type ;
      typedef std::ptrdiff_t                         size_type ;
      typedef contiguous_vector< T, std::ptrdiff_t > base_type ;

    public:
      vector( size_type n )
      : base_type( new T[n], n )
      , owner_( true )
#ifndef NDEBUG
      , exists_( new bool )
#endif
      {
#if !defined(NDEBUG) || defined(GLAS_INITIALIZE_NAN)
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+n, std::numeric_limits<T>::quiet_NaN() ) ;
#elif defined(GLAS_INITIALIZE_ZERO)
        std::fill( this->ptr(), this->ptr()+n, 0.0 ) ;
#endif
#ifndef NDEBUG
        *exists_ = true ;
#endif
      }

      ~vector() {
         if (owner_) {
           delete [] this->ptr() ;
#ifndef NDEBUG
           assert( exists_.use_count()==1 ) ;
#endif
         }
      }

      // Copy reference !!
      vector( vector const& that )
      : base_type( that.ptr(), that.size() )
      , owner_( false )
#ifndef NDEBUG
      , exists_( that.exists_ )
#endif
      {
//        assert( false ) ;
      }
      /*vector( vector const& that )
      : base_type( new T[that.size()], that.size() )
      {
        std::copy( that.ptr(), that.ptr()+that.size(), this->ptr() ) ;
      }*/

      // Move constructor
      vector( vector&& that )
      : base_type( that )
      , owner_( that.owner_ )
#ifndef NDEBUG
      , exists_( new bool(true) )
#endif
      {
        that.owner_ = false ;
#ifdef GLAS2_MEMORY_WARNINGS
        std::cerr << "GLAS2::MATRIX: Move constructor." << std::endl ;
#endif
      }

      // Is a deep copy!!
      vector( std::initializer_list<value_type> that )
      : base_type( new T[that.size()], that.size() )
      , owner_( true )
#ifndef NDEBUG
      , exists_( new bool(true) )
#endif
      {
        std::copy( that.begin(), that.end(), this->ptr() ) ;
      }

    public:
      // Copies of the vector become invalid
      vector& resize( size_type n ) {
        if (owner_) delete [] this->ptr() ;
        this->reset( new T[n], n ) ;
        owner_ = true ;
#ifndef NDEBUG
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+n, std::numeric_limits<T>::quiet_NaN() ) ;
#endif
        return *this ;
      }

      base_type pass_ref() const { return *this ; }

      // Deep copy
      template <typename E>
      vector( copy_of<E> const& that )
      : base_type( new T[that.expression().size()], that.expression().size() )
      , owner_( true )
#ifndef NDEBUG
      , exists_( new bool(true) )
#endif
      {
        *this = that.expression() ;
      }

    public:
      base_type& operator=( vector const& that ) {
#ifndef NDEBUG
        assert( *exists_ ) ;
#endif
        base_type(*this) = that ;
        return *this ;
      }

      base_type& operator=( vector&& that ) = delete ;

    public:
      template <typename E>
      base_type& operator=( E const& that ) {
#ifndef NDEBUG
        assert( *exists_ ) ;
#endif
        return static_cast<base_type&>(*this) = that ;
      }

    private:
      bool owner_ ;
#ifndef NDEBUG
      std::shared_ptr< bool > exists_ ;
#endif
  } ;


  template <typename T>
  struct glas_concept< vector<T> >
  : glas_concept< typename vector<T>::base_type >
  {};

  template <typename T>
  struct pass_reference< vector<T> >
  {
    static const bool pass_by_value = false ;

    typedef typename  vector<T>::base_type type ;

    type operator()( vector<T> const& that ) const {
      return that ;
    }
  };

  template <typename T, typename S2>
  struct vector_selection< vector<T>, S2 >
  : vector_selection< typename vector<T>::base_type, S2 >
  {} ;

} // namespace glas2

#endif
