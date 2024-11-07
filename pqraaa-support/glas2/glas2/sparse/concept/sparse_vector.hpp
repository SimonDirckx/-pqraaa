//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_concept_sparse_vector_hpp
#define glas2_sparse_concept_sparse_vector_hpp

#include <glas2/vector/concept/vector.hpp>

namespace glas2 {

  struct SparseVector
  : Vector
  {
    typedef SparseVector type ;
  } ;

  // concept SparseVector

  // value_type operator()(size_type i)
  // index_type num_nz()
  // size_type size()
  // size_type position(index_type)
  // size_type value(index_type)

} // namespace glas

#endif
