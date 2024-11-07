#ifndef glas2_matrix_algorithm_eye_hpp
#define glas2_matrix_algorithm_eye_hpp

#include <glas2/matrix/type/identity_matrix.hpp>

namespace glas2 {

  inline identity_matrix<double> eye( int m, int n ) { return identity_matrix<double>(m,n) ; }
}

#endif
