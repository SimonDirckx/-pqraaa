#ifndef glas2_tensor_algorithm_tensor_selection_hpp
#define glas2_tensor_algorithm_tensor_selection_hpp

#include <glas2/tensor/concept/dense_tensor.hpp>
#include <cassert>

namespace glas2 {

  template <typename M, typename Selection, typename Mode, typename EnableIf=void>
  struct tensor_selection {
  } ;

}

#endif
