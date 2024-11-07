// Zifan Liu, zifan.liu@cs.kuleuven.be 2015-06-29
//

#ifndef glas2_bsp_matrix_concept_bsp_coord_sparse_matrix_hpp
#define glas2_bsp_matrix_concept_bsp_coord_sparse_matrix_hpp

#include <glas2/bsp/matrix/concept/bsp_matrix.hpp>

namespace glas2 {namespace bsp{

	// BSP distributed sparse matrix
	struct BSPCoordCsrMatrix
	 : BSPMatrix
	{
		typedef BSPCoordCsrMatrix type ;
	} ;
	
	// sparse matrix storage format
	struct COORD_CSR {} ;
	//struct CSC {} ;
	

}}
#endif
