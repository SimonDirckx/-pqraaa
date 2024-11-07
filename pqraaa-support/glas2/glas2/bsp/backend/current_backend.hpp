#ifndef glas2_bsp_backend_current_backend_hpp
#define glas2_bsp_backend_current_backend_hpp

#include <glas2/bsp/backend/default/default_backend.hpp>

namespace glas2 { namespace bsp {

#ifndef GLAS2_BSP_CURRENT_BACKEND
  typedef default_backend current_backend ;
#else
  typedef GLAS2_BSP_CURRENT_BACKEND current_backend ;
#endif

} } // namespace glas2

#endif
