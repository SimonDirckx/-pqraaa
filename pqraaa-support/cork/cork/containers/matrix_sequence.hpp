#ifndef cork_containers_matrix_sequence_hpp
#define cork_containers_matrix_sequence_hpp

#include <glas2/matrix.hpp>

namespace CORK {
  template <typename T, typename S=std::ptrdiff_t, typename O=glas2::column_major>
  using matrix_sequence = glas2::contiguous_matrix<T, S, O>;

  template <typename T>
  class matrix_sequence {

    matrix_sequence( T* ptr,  )

  }

} // namespace CORK

#endif

