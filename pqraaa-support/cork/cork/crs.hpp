//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_crs_hpp
#define cork_crs_hpp


#include <cork/sparse.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <cork/vector.hpp>
#include <type_traits>

namespace CORK {

  template <int IndexBase, typename O, typename T, typename I=int, typename S=I>
  using crs = glas2::crs_structure<T, O, I, S, IndexBase> ;

  template <typename I, typename T>
  using crs0 = glas2::crs_structure<T, glas2::row_major, I, vector_size_type, 0> ;

  template <typename I, typename T>
  using crs1 = glas2::crs_structure<T, glas2::row_major, I, vector_size_type, 1> ;

} // namespace CORK

#endif

