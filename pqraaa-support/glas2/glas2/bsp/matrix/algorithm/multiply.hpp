//Zifan Liu, zifan.liu@gmail.com
#ifndef glas2_bsp_matrix_algorithm_multiply_hpp
#define glas2_bsp_matrix_algorithm_multiply_hpp

#include <glas2/concept/is.hpp>
#include <glas2/matrix/algorithm/multiply.hpp>
#include <glas2/bsp/matrix/concept/bsp_row_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_column_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_csr_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_coord_csr_matrix.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <glas2/bsp/matrix/expression/matrix_vector_multiply.hpp>
//#include <glas2/bsp/matrix/expression/matrix_multiply.hpp>
#include <glas2/bsp/matrix/expression/rmatrix_dvector_multiply.hpp>
#include <glas2/bsp/matrix/expression/cmatrix_vector_multiply.hpp>
#include <glas2/bsp/matrix/expression/csr_mvp.hpp>
#include <glas2/bsp/matrix/expression/csr_tr_mvp.hpp>
#include <glas2/bsp/matrix/expression/coord_csr_mvp.hpp>
#include <glas2/bsp/matrix/expression/coord_csr_mmp.hpp>
#include <type_traits>

//#include <glas2/bsp/backend/current_backend.hpp>
//#include <glas2/bsp/backend/default/matrix/matrix_vector_multiply.hpp>
//#include <glas2/vector/type/contiguous_vector.hpp>

namespace glas2 {


  template <typename M, typename V>
  typename std::enable_if< is<bsp::BSPRowMatrix,M>::value && is<DenseVector,V>::value
                         , bsp::distributed_vector< matrix_vector_multiply<typename M::local_type,V>, typename M::distribution_type >
                         >::type multiply( M const& m, V const& v ) {
    return bsp::distributed_vector< matrix_vector_multiply<typename M::local_type,V>
                                  , typename M::distribution_type
                                  >( multiply( m.local(), v ), m.distribution() ) ;
  }


  // MVP of a column-wise distributed matrix and a distributed vector
/*  template<typename M, typename V>
  typename std::enable_if< is<bsp::BSPRowMatrix, M>::value && is<bsp::BSPVector, V>::value 
			, bsp::distributed_vector<typename glas2::contiguous_vector<typename V::value_type,typename V::size_type>, typename M::distribution_type >
			>::type multiply( M const& m, V const& v) {
	return glas2::bsp::multiply ( bsp::current_backend(), m, v );

			}
*/

   template <typename M, typename V>
   typename std::enable_if< is<bsp::BSPRowMatrix, M>::value && is<bsp::BSPVector, V>::value
			,  bsp::distributed_vector< bsp::rmatrix_dvector_multiply<M,V>, typename M::distribution_type >
		        >::type multiply (M const& m, V const& v) {
	return bsp::distributed_vector< bsp::rmatrix_dvector_multiply<M,V>, typename M::distribution_type >(bsp::rmatrix_dvector_multiply<M,V>(m,v), m.distribution() );
			}	

   template <typename M, typename V>
   typename std::enable_if< is<bsp::BSPColumnMatrix, M>::value && is<DenseVector, V>::value
   			,  bsp::distributed_vector< bsp::cmatrix_vector_multiply<M,V>, typename M::distribution_type >
			>::type multiply (M const& m, V const& v) {
	return bsp::distributed_vector< bsp::cmatrix_vector_multiply<M,V>, typename M::distribution_type >(bsp::cmatrix_vector_multiply<M,V>(m,v), m.distribution() );
			}

   template <typename M, typename V>
   typename std::enable_if< is<bsp::BSPCsrMatrix, M>::value && is<DenseVector, V>::value 
   			, bsp::distributed_vector< bsp::csr_mvp<M,V>, typename M::distribution_type >
			>::type multiply (M const& m, V const& v) {
	return bsp::distributed_vector< bsp::csr_mvp<M,V>, typename M::distribution_type >( bsp::csr_mvp<M,V>(m,v), m.distribution() );
			}

/*   template <typename M, typename V>
   typename std::enable_if<is<bsp::BSPCoordCsrMatrix, M>::value && is<DenseVector, V>::value
   			, bsp::coord_csr_mvp<M,V>
   			>::type multiply (M const& m, V const& v) {
	return bsp::coord_csr_mvp<M,V> (bsp::coord_csr_mvp<M,V>(m,v));
			}
*/
   template <typename M, typename V>
   typename std::enable_if<is<bsp::BSPCoordCsrMatrix, M>::value && is<DenseVector, V>::value
                        , bsp::coord_csr_mvp<M,V>
                        >::type multiply (M const& m, V const& v) {
        return bsp::coord_csr_mvp<M,V>(m,v);
                        }

  template <typename M1, typename M2>
  typename std::enable_if< is<bsp::BSPCoordCsrMatrix,M1>::value && is<DenseMatrix,M2>::value
                         , bsp::coord_csr_mmp<M1,M2>
                         >::type multiply( M1 const& m1, M2 const& m2 ) {
    return bsp::coord_csr_mmp<M1,M2>( m1, m2 ) ;
  }


   template <typename M, typename V>	      
	   typename std::enable_if< is<bsp::BSPCsrMatrix, M>::value && is<DenseVector, V>::value
	   , bsp::distributed_vector< bsp::csr_mvp<M,V>, typename M::distribution_type >
	   >::type transpose_multiply (M const& m, V const& v) {
		   return bsp::distributed_vector< bsp::csr_mvp<M,V>, typename M::distribution_type >( bsp::csr_tr_mvp<M,V>(m,v), m.distribution() );
									                              }
	/*
  template <typename M1, typename M2>
  typename std::enable_if< is<bsp::BSPRowMatrix,M1>::value && is<DenseMatrix,M2>::value
                         , bsp::row_distributed_matrix< matrix_multiply<typename M1::local_type,M2>, typename M1::distribution_type >
                         >::type multiply( M1 const& m1, M2 const& m2 ) {
    return bsp::row_distributed_matrix< matrix_multiply<typename M1::local_type,M2>
                                      , typename M1::distribution_type
                                      >( multiply( m1.local(), m2 ), m1.distribution() ) ;
  }

  template <typename M, typename V>
  typename std::enable_if< is<bsp::BSPColumnMatrix,M>::value && is<bsp::BSPVector,V>::value
                         , bsp::matrix_vector_multiply< M, V >
                         >::type multiply( M const& m, V const& v ) {
    return matrix_vector_multiply< M, V >( m, v ) ;
  }

  template <typename M1, typename M2>
  typename std::enable_if< is<bsp::BSPColumnMatrix,M1>::value && is<bsp::BSPRowMatrix,M2>::value
                         , bsp::matrix_multiply< M1, M2 >
                         >::type multiply( M1 const& m1, M2 const& m2 ) {
    return bsp::matrix_multiply< M1, M2 >( m1, m2 ) ;
  }
*/
} // namespace glas2

#endif
