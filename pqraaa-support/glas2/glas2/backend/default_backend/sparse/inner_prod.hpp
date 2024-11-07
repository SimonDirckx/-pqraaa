//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_sparse_inner_prod_hpp
#define glas2_backend_default_backend_sparse_inner_prod_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/backend/default_backend/matrix/fill.hpp>
#include <glas2/backend/default_backend/sparse/is_coordinate_sparse_matrix_multiply.hpp>
#include <glas2/backend/default_backend/sparse/is_compressed_sparse_matrix_multiply.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/sparse/concept/sparse_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename SV, typename DV>
  typename std::enable_if< is<DenseVector,DV>::value && is< SparseVector, SV >::value
                         , decltype( typename SV::value_type() * typename DV::value_type() )
                         >::type inner_prod( default_backend, SV const& sv, DV const& dv )
  {
    assert( dv.size()==sv.size() );

    decltype( typename SV::value_type() * typename DV::value_type() ) sum = 0.0 ;

    for ( typename SV::size_type j=0; j<sv.num_nz(); ++j ) {
      sum += sv.data()(j) * dv( sv.index()(j) ) ;
    }
    return sum ;
  }

} // namespace glas2

#endif
