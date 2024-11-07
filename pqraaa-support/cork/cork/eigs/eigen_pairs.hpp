//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_eigen_pairs_hpp
#define cork_eigs_eigen_pairs_hpp

#include <cork/vector.hpp>
#include <cork/matrix.hpp>
#include <type_traits>

namespace CORK { namespace eigs {

  template <typename X, typename E>
  class eigen_pairs {
    public:
      typedef typename std::decay<X>::type x_type ;
      typedef typename std::decay<E>::type e_type ;
      typedef typename std::common_type< typename x_type::value_type
                                       , typename e_type::value_type
                                       >::type value_type ;

    public:
      template <typename X1, typename E1>
      eigen_pairs( X1 x_, E1 e_ )
      : x_( x_ )
      , e_( e_ )
      , order_( e_.size() )
      {
        assert( e_.size() == x_.num_columns() ) ;
      }

      eigen_pairs( eigen_pairs const& that )
      : x_( that.x_ )
      , e_( that.e_ )
      , order_( that.order_ )
      {}

      eigen_pairs( eigen_pairs&& that )
      : x_( std::move( that.x_ ) )
      , e_( std::move( that.e_ ) )
      , order_( that.order_ )
      {}

      void size( int ord ) {
        assert( ord<=e_.size() ) ;
        order_ = ord ;
      }

      void operator=( eigen_pairs const& that ) {
        x_ = that.x_ ;
        e_ = that.e_ ;
        order_ = that.order_ ;
      } // move assignment

      void operator=( eigen_pairs&& that ) {
        x_ = std::move( that.x_ ) ;
        e_ = std::move( that.e_ ) ;
        order_ = that.order_ ;
      } // move assignment

    public:
      auto size() const {
        return order_ ;
      }
      auto vectors() const { return CORK::matrix<value_type>(x_.ptr(), x_.num_rows(), order_ ) ; }
      auto values() const { return CORK::vector<value_type>(e_.ptr(), order_ ) ; }

      /*auto vectors() { return CORK::matrix<value_type>(x_.ptr(), x_.num_rows(), order_) ; }
      auto values() { return CORK::vector<value_type>(e_.ptr(), order_ ) ; }*/

    private:
      X x_ ;
      E e_ ;
      int order_ ;
  } ; // class eigen_pairs

} } // namespace CORK::eigs

#endif
