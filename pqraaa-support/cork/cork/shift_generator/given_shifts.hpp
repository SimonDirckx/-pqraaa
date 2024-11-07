//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_shift_generator_given_shifts_hpp
#define cork_shift_generator_given_shifts_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace shift_generator {

  template <typename Shifts>
  class given_shifts_selector
  {
    public:
      typedef typename Shifts::value_type value_type;
      typedef int                         size_type ;

    public:
      explicit given_shifts_selector( Shifts const& shifts )
      : k_( 0 )
      , shifts_( shifts )
      {}
      // delete copy-constructor
      given_shifts_selector(given_shifts_selector const&) = delete;

    public:
      value_type const& next_shift() const {
        int k = k_ ;
        ++k_ ; if (k_>=shifts_.size()) k_ = 0 ;
        return shifts_(k) ;
      } // next_shift()

      value_type const& shift( size_type k ) const {
        return shifts_(k) ; 
      } // shift()
     
      size_type solve_handle() const { return 0 ; }

      template <typename Ritz, typename Resid>
      void update_shifts( Ritz const& ritz_values, Resid const& ) const {}

    private:
      mutable int k_ ;
      Shifts      shifts_ ;
  } ; // class given_shifts_selector

  template <typename Shifts>
  class given_shifts
  {
    public:
      given_shifts( Shifts const& shifts )
      : shifts_( shifts )
      {}

    public:
      template <typename EigenvalueSelector>
      given_shifts_selector<Shifts> operator()(EigenvalueSelector const& eigenvalue_selector) {
        return given_shifts_selector(shifts_);
      }

    private:
      Shifts shifts_ ;
  } ; // class given_shifts

} } // namespace CORK::shift_generator

#endif
