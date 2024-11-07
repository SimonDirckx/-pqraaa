//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_user_defined_eigenvalue_selector_multiple_shifts_hpp
#define cork_user_defined_eigenvalue_selector_multiple_shifts_hpp

#include <cork/shift_generator/eigenvalue_selector_multiple_poles.hpp>

namespace CORK { namespace user_defined {

  class eigenvalue_selector_multiple_shifts
  {
    public:
      explicit eigenvalue_selector_multiple_shifts( int multiplicity )
      : multiplicity_( multiplicity )
      {}

    public:
      int multiplicity() const { return multiplicity_ ; }

      template <typename Geometry>
      shift_generator::eigenvalue_selector_multiple_poles<Geometry> generate_shift_generator( Geometry const& eigenvalue_selector ) {
        return shift_generator::eigenvalue_selector_multiple_poles<Geometry>( eigenvalue_selector, multiplicity_ ) ;
      }

    private:
      int             multiplicity_ ;
  } ; // class eigenvalue_selector_multiple_poles

} } // namespace CORK::user_defined

#endif
