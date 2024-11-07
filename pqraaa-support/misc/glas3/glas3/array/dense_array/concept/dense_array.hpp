//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_concept_dense_array_hpp
#define glas3_array_dense_array_concept_dense_array_hpp

#include <glas3/array/concept/array.hpp>

namespace glas3 {

  struct DenseArray
  : Array
  {
    typedef DenseArray type ;
  } ;

  struct DenseMatrix
  : Matrix, DenseArray
  {
    typedef DenseMatrix type ;
  } ;

  struct DenseVector
  : Vector, DenseArray
  {
    typedef DenseVector type ;
  } ;

  struct DenseScalar
  : Scalar, DenseArray
  {
    typedef DenseScalar type ;
  } ;

} // namespace glas3

#endif
