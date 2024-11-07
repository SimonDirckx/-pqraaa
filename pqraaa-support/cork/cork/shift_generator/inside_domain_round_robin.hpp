//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_shift_generator_inside_domain_round_robin_hpp
#define cork_shift_generator_inside_domain_round_robin_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#include <iostream>

namespace CORK { namespace shift_generator {

  template <typename EigenvalueSelector>
  class inside_domain_round_robin_selector
  {
    public:
      typedef typename EigenvalueSelector::value_type     value_type ;
      typedef int                                         size_type ;

    public:
      explicit inside_domain_round_robin_selector(EigenvalueSelector const& eig_selector)
      : shifts_(eig_selector.domain().discretize(eig_selector.multiplicity()))
      , k_(0)
      { }

      // delete copy-constructor
      inside_domain_round_robin_selector(inside_domain_round_robin_selector&) = delete;

    public:
      value_type const& next_shift() const {
        ++k_ ;
        return shift(k_) ;
      } // next_shift()

      value_type const& shift( size_type k ) const {
        return shifts_( k % shifts_.size() ) ;
      } // shift()

      size_type solve_handle() const { return 0 ; }

      template <typename Ritz, typename Resid>
      void update_shifts( Ritz const& ritz_values, Resid const& ) const {}

    private:
      glas2::shared_vector<value_type>   shifts_            ;
      mutable int                        k_                 ;
  } ; // class inside_domain_round_robin_selector

  class inside_domain_round_robin
  {
    public:
      template <typename EigenvalueSelector>
      inside_domain_round_robin_selector<EigenvalueSelector> operator()(EigenvalueSelector const& eigenvalue_selector) {
        return inside_domain_round_robin_selector<EigenvalueSelector>(eigenvalue_selector);
      }
  } ; // class inside_domain_round_robin

} } // namespace CORK::shift_generator

#endif
