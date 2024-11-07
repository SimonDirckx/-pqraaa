//  (C) Copyright Karl Meerbergen 2011. 
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_iterative_krylov_tolerance_hpp
#define glas2_iterative_krylov_tolerance_hpp

#include <glas2/iterative/krylov/options.hpp>

namespace glas2 { namespace iterative {

  template <typename T>
  double tolerance ( options const& opt, T const& norm_b ) {
    return opt.absolute_tolerance_ + opt.relative_tolerance_ * norm_b ;
  }

} }

#endif
