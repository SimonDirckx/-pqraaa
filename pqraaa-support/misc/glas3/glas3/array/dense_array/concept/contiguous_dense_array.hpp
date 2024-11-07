//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_concept_contiguous_dense_array_hpp
#define glas3_array_dense_array_concept_contiguous_dense_array_hpp

#include <glas3/array/dense_array/concept/dense_array.hpp>

namespace glas3 {

  struct ContiguousDenseArray
  : DenseArray
  {
    typedef ContiguousDenseArray type ;
  } ;

  struct ContiguousDenseMatrix
  : DenseMatrix, ContiguousDenseArray
  {
    typedef ContiguousDenseMatrix type ;
  } ;

  struct ContiguousDenseVector
  : DenseVector, ContiguousDenseArray
  {
    typedef ContiguousDenseVector type ;
  } ;

  struct ContiguousDenseScalar
  : DenseScalar, ContiguousDenseArray
  {
    typedef ContiguousDenseScalar type ;
  } ;

} // namespace glas3

#endif
