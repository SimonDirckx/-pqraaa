//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_scalar_container_shared_scalar_hpp
#define glas2_scalar_container_shared_scalar_hpp

#include <type_traits>
#include <memory>

namespace glas2 {

  template <typename T>
  class shared_scalar
  {
    public:
      typedef T                         value_type ;

    public:
      shared_scalar( T const& value=T() )
      : shared_( new T )
      {
        *shared_ = value ;
      }

      // Copy reference !!
      shared_scalar( shared_scalar const& that )
      : shared_( that.shared_ )
      {}

    public:
      shared_scalar& operator=( shared_scalar const& that ) {
        *shared_ = *(that.shared_) ;
        return *this ;
      }

      shared_scalar& operator=( shared_scalar&& that ) = delete ;

      // Deep copy
      shared_scalar copy() const {
        shared_scalar that ;
        that = *this ;
        return that ;
      }

    public:
      T& value() { return *shared_ ; }
      T const& value() const { return *shared_ ; }

      operator T() const { return *shared_ ; }
      operator T&() { return *shared_ ; }

    public:
      template <typename E>
      shared_scalar& operator=( E const& that ) {
        *shared_ = that ;
        return *this ;
      }

    private:
      std::shared_ptr<value_type> shared_ ;
  } ;


} // namespace glas2

namespace std {

  template <typename T>
  struct is_integral< glas2::shared_scalar<T> >
  : is_integral<T>
  {} ;

  template <typename T>
  struct is_arithmetic< glas2::shared_scalar<T> >
  : is_arithmetic<T>
  {} ;

}

#endif
