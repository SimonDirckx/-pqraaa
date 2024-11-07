//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coo_hpp
#define cork_coo_hpp


#include <cork/sparse.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <cork/vector.hpp>
#include <type_traits>

namespace CORK {

  template <int IndexBase, typename I, typename T>
  using coo = glas2::coo_structure<T, I, vector_size_type, IndexBase> ;

  template <typename I, typename T>
  using coo0 = glas2::coo_structure<T, I, vector_size_type, 0> ;

  template <typename I, typename T>
  using coo1 = glas2::coo_structure<T, I, vector_size_type, 1> ;

} // namespace CORK

#endif

