//  (C) Copyright Karl Meerbergen, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_concept_hpp
#define cork_matrix_valued_function_concept_hpp

#include <cork/concept/config.hpp>
#include <cork/concept/vector_of.hpp>

#ifdef CORK_USE_CONCEPTS

#include <cork/vector.hpp>
#include <type_traits>
#include <concepts>
#include <cmath>

namespace CORK {

  template <typename MVF>
  concept MatrixValuedFunction = requires( MVF& mvf
                                         , typename MVF::size_type  n
                                         , typename MVF::size_type i
                                         , CORK::vector<typename MVF::value_type> x
                                         , CORK::vector<typename MVF::value_type> y
                                         , typename MVF::value_type& shift
                                         , typename MVF::value_type alpha
                                         , bool is_new_shift
                                         , typename MVF::size_type clone ) {
    typename MVF::value_type ;
    typename MVF::size_type ;
    {d.matvec(i,alpha,x,y)} -> std::same_as<void> ;
    {d.solve(shift, clone, x, is_new_shift)} -> std::same_as<void> ;
    {d.size()} -> std::same_as< size_type > ;
    {d.num_terms()} -> std::same_as< size_type > ;
  } ;

} // namespace CORK

#endif

#endif

