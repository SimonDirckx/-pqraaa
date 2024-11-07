#ifndef glas2_backend_current_backend_hpp
#define glas2_backend_current_backend_hpp

#include <glas2/backend/default_backend/default_backend.hpp>

namespace glas2 {

#ifndef GLAS2_CURRENT_BACKEND
  typedef default_backend current_backend ;
#else
  typedef GLAS2_CURRENT_BACKEND current_backend ;
#endif

} // namespace glas2

#endif
