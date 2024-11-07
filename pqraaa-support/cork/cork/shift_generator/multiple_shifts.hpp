//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_shift_generator_multiple_shifts_hpp
#define cork_shift_generator_multiple_shifts_hpp

#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>
#include <cmath>

namespace CORK { namespace shift_generator {

  template <typename Shifts>
  class multiple_shifts
  {
    public:
      typedef typename std::decay<Shifts>::type::value_type value_type ;
      typedef typename std::decay<Shifts>::type::size_type  size_type ;

    public:
      explicit multiple_shifts( Shifts shifts, int multiplicity )
      : shifts_( shifts )
      , multiplicity_( multiplicity )
      , k_(0)
      {
        assert( shifts.size()>0 ) ;
      }

    public:
      value_type next_shift() const {
        ++k_ ;
        value_type shift = shifts_( ((k_-1)/multiplicity_) % shifts_.size() ) ;
        for ( ; std::isnan(glas2::real(shift)); ++k_ ) {
          shift = shifts_( ((k_-1)/multiplicity_) % shifts_.size() ) ;
        }
        return shift ;
      }

      size_type solve_handle() const {
        return 0 ;
      }

      value_type const& shift( size_type k ) const {
        return shifts_((k/multiplicity_) % shifts_.size() ) ;
      }

      void modify_shift( value_type const& shift ) {
        shifts_( ((k_-1)/multiplicity_) % shifts_.size() ) = shift ;
      }

    public:
      template <typename Ritz>
      void update_shifts( Ritz const& ritz_values ) const {}

      template <typename Ritz, typename Resid>
      void update_shifts( Ritz const& ritz_values, Resid const& resid ) const {}

    private:
      Shifts       shifts_ ;
      int         multiplicity_ ;
      mutable int k_ ;
  } ; // class multiple_shifts

} } // namespace CORK::shift_generator

#endif
