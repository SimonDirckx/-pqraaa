#ifndef glas2_sparse_algorithm_speye_hpp
#define glas2_sparse_algorithm_speye_hpp

#include <glas2/sparse/type/sparse_identity_matrix.hpp>

namespace glas2 {

  inline sparse_identity_matrix<double> speye( int m, int n ) { return sparse_identity_matrix<double>(m,n) ; }
}

#endif
