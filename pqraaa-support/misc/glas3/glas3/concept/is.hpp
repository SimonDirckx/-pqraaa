//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_concept_is_hpp
#define glas3_concept_is_hpp

#include <glas3/concept/concept.hpp>
#include <type_traits>

namespace glas3 {

//struct Vector ;

template <typename Concept, typename Class, typename EnableIf=void>
struct is
//: boost::is_base_of< traits<Class>::concept, Concept >
: std::is_base_of< Concept, typename concept<Class>::type >
{} ;

template <typename Concept, typename Class>
struct is<Concept, Class&>
: is<Concept, Class>
{} ;

template <typename Concept, typename Class>
struct is<Concept, Class const>
: is<Concept, Class>
{} ;

//template<class T, class...>
//struct are_same : std::true_type
//{};
//
//template<class T, class U, class... TT>
//struct are_same<T, U, TT...>
//: std::integral_constant<bool, std::is_same<T,U>{} && are_same<T, TT...>{}>
//{};

template<typename... I>
struct are_integral
: std::true_type
{};

template<typename I, typename... J>
struct are_integral<I, J...>
: std::integral_constant<bool, std::is_integral<I>::value && are_integral<J...>::value>
{};

} // namespace glas3

#endif
