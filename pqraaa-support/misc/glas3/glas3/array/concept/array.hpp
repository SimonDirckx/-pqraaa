//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_concept_array_hpp
#define glas3_array_concept_array_hpp

namespace glas3 {

  struct Array
  {
    typedef Array type ;
  } ;

  // concept Array:
  //
  // value_type
  // shape_type
  // size_type
  // ndims_type
  // shape()
  // size()
  // [i] -> i-th entry

  struct Matrix
  : Array
  {
    typedef Matrix type ;
  } ;

  struct Vector
  : Array
  {
    typedef Vector type ;
  } ;

  struct Scalar
  : Array
  {
    typedef Scalar type ;
  } ;
//
//  struct EmptyArray
//  : Array
//  {
//    typedef EmptyArray type ;
//  } ;


} // namespace glas3

#endif
