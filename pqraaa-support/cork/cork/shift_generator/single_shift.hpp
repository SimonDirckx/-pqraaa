//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_shift_generator_single_shift_hpp
#define cork_shift_generator_single_shift_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace shift_generator {

  template <typename ShiftValueType>
  class single_shift_selector
  {
    public:
      typedef ShiftValueType     value_type;
      typedef int                size_type ;

    public:
      explicit single_shift_selector( value_type const& shift )
      : shift_( shift )
      {}

      // delete copy-constructor
      single_shift_selector(single_shift_selector const&) = delete;

    public:
      value_type const& next_shift() const {
        return shift_ ;
      } // next_shift()

      value_type const& shift( size_type k ) const {
        return shift_ ; 
      } // shift()
     
      size_type solve_handle() const { return 0 ; }

      template <typename Ritz, typename Resid>
      void update_shifts( Ritz const& ritz_values, Resid const& ) const {}

    public:
      template <typename Ritz>
      void update_shifts( Ritz const& ritz_values ) const {}

      void modify_shift( value_type const& shift ) {
        shift_ = shift ;
      }

    private:
      value_type shift_ ;
  } ; // class single_shift_selector

  class single_shift
  {
    public:
      template <typename EigenvalueSelector>
      single_shift_selector<EigenvalueSelector> operator()(EigenvalueSelector const& eigenvalue_selector) {
        return single_shift_selector<EigenvalueSelector>(eigenvalue_selector);
      }
  } ; // class single_shift

} } // namespace CORK::shift_generator

#endif
