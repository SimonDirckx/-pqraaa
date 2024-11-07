//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_barycentric_rational_hpp
#define cork_barycentric_rational_hpp

#include <cork/vector.hpp>
#include <cork/basis/barycentric_rational.hpp>
#include <cork/basis4cork/barycentric_rational.hpp>
#include <cork/matrix_iterator/barycentric_rational.hpp>

namespace CORK {

  template <typename T>
  using barycentric_rational = basis::barycentric_rational< CORK::vector<T>, CORK::vector<T> > ;

} // namespace CORK

#endif
