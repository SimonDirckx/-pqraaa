//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_concept_stl_sequence_container_hpp
#define cork_concept_stl_sequence_container_hpp

#include <type_traits>
#include <cork/concept/config.hpp>

#ifdef CORK_USE_CONCEPTS
#include <concepts>
#endif

namespace CORK {

#ifdef CORK_USE_CONCEPTS
  template <typename T>
  concept STLForwardIterator = requires ( T it, T it2 ) {
    {it==it2} -> std::same_as<bool> ;
    {it!=it2} -> std::same_as<bool> ;
    {++it} ;
    {it++} ;
    {T(it)} ;
    {it=it2} -> std::same_as<T&> ;
    {*it} ;
  } ;

  template <typename F>
  concept STLSequenceContainer = requires ( F const c_sequence, F sequence ) {
    typename F::value_type ;
    typename F::size_type ;
    typename F::iterator ;
    typename F::const_iterator ;

    requires STLForwardIterator< typename F::iterator > ;
    requires STLForwardIterator< typename F::const_iterator > ;

    {c_sequence.size()} -> std::same_as< typename F::size_type > ;
    {c_sequence.begin()} -> std::same_as< typename F::const_iterator > ;
    {c_sequence.end()} -> std::same_as< typename F::const_iterator > ;

    {sequence.begin()} -> std::same_as< typename F::iterator > ;
    {sequence.end()} -> std::same_as< typename F::iterator > ;
  } ;
#endif

} // namespace CORK

#endif

