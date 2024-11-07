//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_ode_filter_unstable_hpp
#define cork_matrix_ode_filter_unstable_hpp

#include <cork/options/filter_poles.hpp>
#include <cork/vector.hpp>

namespace CORK { namespace ode {

  namespace detail {
    template <typename T>
    int filter_unstable( CORK::vector<T> eig_val ) {
      int n_keep = 0 ;
      for (int i=0; i<eig_val.size(); ++i) {
        if (eig_val(i).real()<=0.0) {
          eig_val(n_keep) = eig_val(i) ;
          ++n_keep ;
        }
      }
      return n_keep ;
    }
  } // namespace detail

  template <typename T>
  class filter_unstable
  : public options::filter_poles<T> {
    public:
      filter_unstable()
      : options::filter_poles<T>( detail::filter_unstable<T> )
      {}
  } ;

} } // namespace CORK::ode

#endif
