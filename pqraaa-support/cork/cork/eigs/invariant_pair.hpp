//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_invariant_pair_hpp
#define cork_eigs_invariant_pair_hpp

#include <type_traits>

namespace CORK { namespace eigs {

  template <typename X, typename S, typename T>
  class invariant_pair {
    public:
      typedef typename std::decay<X>::type x_type ;
      typedef typename std::decay<S>::type s_type ;
      typedef typename std::decay<T>::type t_type ;
      typedef typename std::common_type< typename x_type::value_type
                                       , typename std::common_type< typename s_type::value_type, typename t_type::value_type>::type
                                       >::type value_type ;

    public:
      template <typename XX, typename SS, typename TT>
      invariant_pair( XX x_, SS s_, TT t_, int order )
      : x_( x_ )
      , s_( s_ )
      , t_( t_ )
      , order_( order )
      {
        assert( s_.num_columns() == t_.num_columns() ) ;
        assert( order<=s_.num_columns() ) ;
      }

      invariant_pair( invariant_pair const& that )
      : x_( that.x_ )
      , s_( that.s_ )
      , t_( that.t_ )
      , order_( that.order_ )
      {}

      void order( int ord ) {
        assert( ord<=s_.num_columns() ) ;
        order_ = ord ;
      }

      void operator=( invariant_pair const& that ) {
        x_ = that.x_ ;
        s_ = that.s_ ;
        t_ = that.t_ ;
        order_ = that.order_ ;
      } // move assignment

      void operator=( invariant_pair&& that ) {
        x_ = std::move( that.x_ ) ;
        s_ = std::move( that.s_ ) ;
        t_ = std::move( that.t_ ) ;
        order_ = that.order_ ;
      } // move assignment

    public:
      auto vectors() const { return x_( glas2::all(), glas2::range(0,order_) ) ; }
      auto matrix_S() const { return s_( glas2::range(0,order_), glas2::range(0,order_) ) ; }
      auto matrix_T() const { return t_( glas2::range(0,order_), glas2::range(0,order_) ) ; }

    private:
      X x_ ;
      S s_ ;
      T t_ ;
      int order_ ;
  } ;

} } // namespace CORK::eigs

#endif
