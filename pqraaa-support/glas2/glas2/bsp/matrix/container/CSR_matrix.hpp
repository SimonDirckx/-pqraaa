// Zifan Liu, zifan.liu@cs.kuleuven.be


#ifndef glas2_bsp_matrix_container_sparse_matrix_hpp
#define glas2_bsp_matrix_container_sparse_matrix_hpp

#include <glas2/vector/type/contiguous_vector.hpp>

#include <glas2/bsp/matrix/concept/bsp_csr_matrix.hpp>
#include <type_traits>
#include <cassert>
#include <iostream>


#endif //glas2_bsp_matrix_container_sparse_matrix_hpp

namespace glas2 { namespace bsp {

	template< typename Distribution, typename T, typename F=CSR, typename S=std::ptrdiff_t>
	class BSP_CSR_matrix {
	public :
	typedef contiguous_vector<T,S>	base_type;
	typedef typename base_type::value_type 	value_type;
	typedef typename base_type::size_type	size_type;
	typedef contiguous_vector<size_type, S>	index_array_type;
	typedef Distribution distribution_type;

	private : 
	// The sparse matrix is stored in CSR format
	
	/** Array containing the values of the matrix. */
 	base_type val;

	/** Array containing the column index of nonzeros */
	index_array_type col_ind;

	/** Array conainint the poniter to the begining of each row */
	index_array_type row_ptr;

	Distribution d;


	public :

	/**
	 ** Base constructor.
	 **/
	BSP_CSR_matrix(size_type Nnz, Distribution d_)
	: val (new value_type[Nnz], Nnz ),
	col_ind ( new size_type[Nnz],Nnz),
	row_ptr ( new size_type[d_.getSize()+1], d_.getSize()+1 ),
	d(d_)
	{	
		// Initialization
		for (size_type i = 0; i < d.getSize(); ++i)
			row_ptr(i) = -1; // If the row is not yet entered, we mark its start point as -1
		row_ptr(d.getSize()) = 0; // we calculate the nnz value for this processor during the fill-in step 
	}
	
	/**
	 ** Initialize a CSR matrix from coordinate format file
	 **/
	void operator() ( size_type i, size_type j, value_type v) {
		size_type siz = row_ptr(d.getSize());
		//if (bsp_pid() == 2)
		//	std::cout<<siz<<std::endl;
		val(siz) = v;
		col_ind(siz) = j;
		size_type ind = i - d.getMinLocalGid();
		if (row_ptr(ind) == -1)
			row_ptr(ind) = siz;
		row_ptr(d.getSize()) ++;
	}

	Distribution const& distribution() const { return d ; }

	value_type getVal (size_type i) const {
		return val(i);
	}

	size_type getCol_ind (size_type i) const {
		return col_ind(i);
	}

	size_type getRow_ptr (size_type i) const {
		return row_ptr(i);
	}

	size_type getNNZ () const {
		return row_ptr(d.getSize());
	}

	};
}}

namespace glas2 {
template <typename D, typename T, typename F>
struct concept< bsp::BSP_CSR_matrix<D,T,F> >
: bsp::BSPCsrMatrix
{};

}

