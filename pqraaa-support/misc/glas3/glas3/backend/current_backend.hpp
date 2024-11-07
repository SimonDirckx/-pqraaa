//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_backend_current_backend_hpp
#define glas3_backend_current_backend_hpp

#include <glas3/backend/default_backend/default_backend.hpp>

namespace glas3 {

#ifndef GLAS3_CURRENT_BACKEND
  typedef default_backend current_backend ;
#else
  typedef GLAS3_CURRENT_BACKEND current_backend ;
#endif

} // namespace glas3

#endif
