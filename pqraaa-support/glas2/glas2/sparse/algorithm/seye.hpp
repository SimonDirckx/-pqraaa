#ifndef glas2_sparse_algorithm_seye_hpp
#define glas2_sparse_algorithm_seye_hpp

#include <glas2/sparse/type/sparse_identity_matrix.hpp>

namespace glas2 {

  inline sparse_identity_matrix<double> seye( int m, int n ) { return sparse_identity_matrix<double>(m,n) ; }
}

#endif
