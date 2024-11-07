//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_container_matrix_hpp
#define glas2_matrix_container_matrix_hpp

#include <glas2/type/pass_reference.hpp>
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
  class matrix
  : public contiguous_matrix< T, std::ptrdiff_t, O >
  {
    public:
      typedef contiguous_matrix< T, std::ptrdiff_t, O > base_type ;
      typedef T              value_type ;
      typedef std::ptrdiff_t size_type ;

    public:
      matrix( size_type m, size_type n )
      : contiguous_matrix< T, std::ptrdiff_t, O >( new T[m*n], m, n )
      , owner_( true )
#ifndef NDEBUG
      , valid_( true )
#endif
#ifndef NDEBUG
      , exists_( new bool(true) )
#endif
      {
#if !defined(NDEBUG) || defined(GLAS_INITIALIZE_NAN)
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+m*n, std::numeric_limits<T>::quiet_NaN() ) ;
#elif defined(GLAS_INITIALIZE_ZERO)
        std::fill( this->ptr(), this->ptr()+m*n, 0.0 ) ;
#endif
      }

      // Shallow Copy !!
      matrix( matrix const& that )
      : contiguous_matrix< T, std::ptrdiff_t, O >( that.ptr(), that.num_rows(), that.num_columns() )
      , owner_( false )
#ifndef NDEBUG
      , valid_( true )
#endif
#ifndef NDEBUG
      , exists_( that.exists_ )
#endif
      {}

      // Move constructor
      matrix( matrix&& that )
      : contiguous_matrix< T, std::ptrdiff_t, O >( that.ptr(), that.num_rows(), that.num_columns() )
      , owner_( that.owner_ )
#ifndef NDEBUG
      , valid_( true )
#endif
#ifndef NDEBUG
      , exists_( new bool(true) )
#endif
      {
        that.owner_ = false ;
#if defined(GLAS2_MEMORY_WARNINGS) && !defined(NDEBUG)
        std::cerr << "GLAS2::MATRIX: Move constructor." << std::endl ;
        that.valid_ = false ;
#endif
      }

      ~matrix() {
        if (owner_) {
           delete [] this->ptr() ;
#ifndef NDEBUG
           assert( exists_.use_count()==1 ) ;
#endif
        }
      }

    public:
      // obsolete.
      base_type const& pass_ref() const{ 
        assert( valid_ ) ;
        return *this; 
      }

      // Copies of the matrix become invalid
      matrix& resize( size_type n, size_type m ) {
        if (owner_) delete [] this->ptr() ;
        this->reset( new T[n*m], n, m ) ;
        owner_ = true ;
#ifndef NDEBUG
        if (std::numeric_limits<T>::has_quiet_NaN)
          std::fill( this->ptr(), this->ptr()+n*m, std::numeric_limits<T>::quiet_NaN() ) ;
#endif
        return *this ;
      }

      friend void swap( matrix& m1, matrix& m2 ) {
        swap( static_cast<base_type&>(m1), static_cast<base_type&>(m2) ) ;
        bool own = m1.owner_ ; m1.owner_ = m2.owner_ ; m2.owner_ = own ;
#ifndef NDEBUG
        bool val = m1.valid_ ; m1.valid_ = m2.valid_ ; m2.valid_ = val ;
        m1.exists_.swap( m2.exists_ ) ;
#endif
      }

    public:
      base_type& operator=( matrix const& that ) {
        assert( valid_ ) ;
        static_cast<base_type&>(*this) = that ;
        return *this ;
      }

      template <typename E>
      base_type& operator=( E const& that ) {
        assert( valid_ ) ;
        static_cast<base_type&>(*this) = that ;
        return *this ;
      }

      /*template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< base_type, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) {
        return matrix_selection< base_type, I1, I2 >::apply( *this, s1, s2 ) ;
      }*/

    private:
      bool owner_ ;
#ifndef NDEBUG
      bool valid_ ;
#endif
#ifndef NDEBUG
      std::shared_ptr< bool > exists_ ;
#endif
  } ;


  template <typename T, typename O>
  struct glas_concept< matrix<T,O> >
  : glas_concept< typename matrix<T,O>::base_type >
  {};

  template <typename Tag, typename T, typename O>
  struct matrix_transform< Tag, matrix<T,O> >
  : matrix_transform< Tag, typename matrix<T,O>::base_type >
  {} ;

  template <typename T, typename O, typename S1, typename S2>
  struct matrix_selection< matrix<T,O>, S1, S2 >
  : matrix_selection< typename matrix<T,O>::base_type, S1, S2 >
  {} ;

  template <typename T, typename O>
  struct pass_reference< matrix<T,O> >
  {
    static const bool pass_by_value = false ;

    typedef typename  matrix<T,O>::base_type type ;

    type operator()( matrix<T,O> const& that ) const {
      return that ;
    }
  };

} // namespace glas2

#endif
