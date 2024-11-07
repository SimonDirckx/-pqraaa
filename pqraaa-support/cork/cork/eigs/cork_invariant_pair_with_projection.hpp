//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_cork_invariant_pair_hpp
#define cork_eigs_cork_invariant_pair_hpp

//#include <cork/eigs/detail/sort_ritz_values.hpp>
//#include <cork/eigs/detail/implicit_restart.hpp>
#include <cork/eigs/explicit_projection.hpp>
#include <cork/eigs/invariant_pair.hpp>
#include <cork/eigs/invariant_pair_for_projection.hpp>
//#include <cork/eigs/info.hpp>
//#include <cork/options/value_of.hpp>
//#include <cork/options/explicit_projection.hpp>

namespace CORK { namespace eigs {

 template <typename CorkResult, typename Options>
 auto cork_invariant_pair( CorkResult const& cork_result, Options const& options ) {
   typedef typename CorkResult::quadruple_type     quadruple_type ;
   typedef typename quadruple_type::value_type     value_type ;
   typedef decltype(std::abs(value_type()))        real_value_type ;

   auto& information = cork_result.information ;

   auto inv_pair = make_invariant_pair_for_projection( cork_result.quadruple, information.number_converged ) ;

   if (options::value_of<options::explicit_projection>( options )) {
     inv_pair = explicit_projection( cork_result.linearization.matrix_polynomial(), cork_result.quadruple, cork_result.eigenvalue_selector, options,  information ) ;
   }

   return inv_pair ;
 } // cork_invariant_pair()

} } // namespace CORK::eigs


#endif
