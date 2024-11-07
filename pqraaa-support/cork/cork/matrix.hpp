//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_hpp
#define cork_matrix_hpp

#include <glas2/matrix.hpp>

namespace CORK {
  template <typename T, typename S=std::ptrdiff_t, typename O=glas2::column_major>
  using matrix = glas2::contiguous_matrix<T, S, O>;

} // namespace CORK

#endif

