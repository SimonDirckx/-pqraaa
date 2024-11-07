//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_rational_newton_hpp
#define cork_rational_newton_hpp

#include <cork/vector.hpp>
#include <cork/basis/rational_newton.hpp>
#include <cork/basis4cork/rational_newton.hpp>
#include <cork/matrix_iterator/rational_newton.hpp>

namespace CORK {

  template <typename T>
  using rational_newton = basis::rational_newton< CORK::vector<T> > ;

} // namespace CORK::basis

#endif
