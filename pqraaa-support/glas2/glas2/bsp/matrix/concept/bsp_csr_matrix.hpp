// Zifan Liu, zifan.liu@cs.kuleuven.be 2015-05-19
//

#ifndef glas2_bsp_matrix_concept_bsp_sparse_matrix_hpp
#define glas2_bsp_matrix_concept_bsp_sparse_matrix_hpp

#include <glas2/bsp/matrix/concept/bsp_matrix.hpp>

namespace glas2 {namespace bsp{

	// BSP distributed sparse matrix
	struct BSPCsrMatrix
	 : BSPMatrix
	{
		typedef BSPCsrMatrix type ;
	} ;
	
	// sparse matrix storage format
	struct CSR {} ;
	//struct CSC {} ;
	

}}
#endif
