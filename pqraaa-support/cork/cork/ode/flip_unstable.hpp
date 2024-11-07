//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_ode_flip_unstable_hpp
#define cork_matrix_ode_flip_unstable_hpp

#include <cork/options/filter_poles.hpp>
#include <cork/vector.hpp>

namespace CORK { namespace ode {

  namespace detail {
    template <typename T>
    int flip_unstable( CORK::vector<T> eig_val ) {
      for (int i=0; i<eig_val.size(); ++i) {
        if (eig_val(i).real()>0.0) {
          eig_val(i) = -std::conj(eig_val(i)) ;
        }
      }
      return eig_val.size() ;
    }
  } // namespace detail

  template <typename T>
  class flip_unstable
  : public options::filter_poles<T> {
    public:
      flip_unstable()
      : options::filter_poles<T>( detail::flip_unstable<T> )
      {}
  } ;

} } // namespace CORK::ode

#endif
