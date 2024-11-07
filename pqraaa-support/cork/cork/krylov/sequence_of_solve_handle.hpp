//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_sequence_of_solve_handle_hpp
#define cork_krylov_sequence_of_solve_handle_hpp

#include <cork/linearization/cork_linearization.hpp>
#include <vector>

namespace CORK { namespace krylov {

  template <typename ValueType, typename Linearization>
  class sequence_of_solve_handle
  {
    public:
      typedef Linearization                                                        linearization_type ;
      typedef typename linearization_type::template solve_handle_type< ValueType > solve_handle_type ;
      typedef typename solve_handle_type::value_type                               value_type ;

    private:
      typedef std::vector< solve_handle_type >  sequence_type ;
      typedef typename sequence_type::size_type size_type ;

    public:
      sequence_of_solve_handle( linearization_type& cork_linearization )
      : linearization_( cork_linearization )
      , number_(0)
      {
        sequence_.push_back( solve_handle_type(linearization_.template solve_handle<ValueType>()) ) ;
        initialized_.push_back( false ) ;
      }

    public:
      auto grade() const { return linearization_.grade() ; }
      auto const& linear_solver() const {
        assert( initialized_[number_] ) ;
        return sequence_[number_].linear_solver() ;
      }

      decltype(auto) basis_handle() const {
        assert( initialized_[number_] ) ;
        return sequence_[number_].basis_handle() ;
      }

    public:
      void select( size_type number ) {
        assert( number>=0 && number<=sequence_.size() ) ;
        if (number==sequence_.size()) {
          sequence_.push_back( linearization_.template solve_handle<ValueType>() ) ;
          initialized_.push_back( false ) ;
        }
        number_ = number ;
      }

      void shift( value_type s ) {
        if (!initialized_[number_] || s!=shift()) {
          sequence_[number_].shift( s ) ;
        }
        initialized_[number_] = true ;
      }

      value_type shift() const {
        assert( initialized_[number_] ) ;
        return sequence_[number_].shift() ;
      }

      template <typename QQ, typename UUIn, typename WOut, typename UUOut, typename Backend>
      void solve_upper( QQ const& Q, UUIn const& U_in, WOut w_out, UUOut U_out, Backend const& backend ) {
        assert( initialized_[number_] ) ;
        sequence_[number_].solve_upper( Q, U_in, w_out, U_out, backend ) ;
      }

      template <typename UUOut>
      void solve_lower( UUOut U_out ) {
        assert( initialized_[number_] ) ;
        sequence_[number_].solve_lower( U_out ) ;
      }

    private:
      linearization_type& linearization_ ;
      sequence_type       sequence_ ;
      size_type           number_ ;
      std::vector<bool>   initialized_ ;
  } ;


} } // namespace CORK::krylov

#endif
