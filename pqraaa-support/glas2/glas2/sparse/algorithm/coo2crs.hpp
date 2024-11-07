#ifndef glas2_sparse_algorithm_coo2crs_hpp
#define glas2_sparse_algorithm_coo2crs_hpp

namespace glas2 {

  template <typename COO, typename CSR>
  CRS& coo2crs( COO& coo, CRS& crs ) {
    crs.reset( coo.num_rows(), coo.num_columns(), coo.nnz() ) ;
    return crs ;
  }

}

#endif
