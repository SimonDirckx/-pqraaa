// Zifan Liu, zifan.liu@cs.kuleuven.be


#ifndef glas2_bsp_matrix_container_coord_sparse_matrix_hpp
#define glas2_bsp_matrix_container_coord_sparse_matrix_hpp

#include <glas2/vector/type/contiguous_vector.hpp>

#include <glas2/bsp/matrix/concept/bsp_coord_csr_matrix.hpp>
#include <type_traits>
#include <cassert>
#include <iostream>


#endif //glas2_bsp_matrix_container_coord_sparse_matrix_hpp

// The difference from the BSP_CSR_matrix is that this class do not need val array, since each of its entry will be zero.

namespace glas2 { namespace bsp {

	//template< typename Distribution, typename T, typename F=CSR, typename S=std::ptrdiff_t>
	template< typename Distribution, typename F=COORD_CSR, typename S=std::ptrdiff_t>
	class BSP_COORD_CSR_matrix {
	public :
	// typedef contiguous_vector<T,S>	base_type;
	 typedef unsigned int 	value_type;
	typedef S	size_type;
	typedef contiguous_vector<size_type, S>	index_array_type;
	typedef Distribution distribution_type;

	private : 
	// The sparse matrix is stored in COORD CSR format
	
	/** Array containing the values of the matrix. */
 	// base_type val;

	/** Array containing the column index of nonzeros */
	index_array_type col_ind;

	/** Array conainint the poniter to the begining of each row */
	index_array_type row_ptr;

	Distribution const& d;

	size_type const& nbrows_;
	size_type const& nbcols_;


	public :

	/**
	 ** Base constructor.
	 **/
	BSP_COORD_CSR_matrix(size_type const& Nnz, Distribution const& d_ , size_type const& nbrows, size_type const& nbcols)
	: col_ind ( new size_type[Nnz],Nnz),
	row_ptr ( new size_type[d_.getSize()+1], d_.getSize()+1 ),
	d(d_),
	nbrows_(nbrows),
	nbcols_(nbcols)
	{	
		// Initialization
		for (size_type i = 0; i < d.getSize(); ++i)
			row_ptr(i) = -1; // If the row is not yet entered, we mark its start point as -1
		row_ptr(d.getSize()) = 0; // we calculate the nnz value for this processor during the fill-in step 
	}
	
	/**
	 ** Initialize a COORD CSR matrix from coordinate format file
	 **/
	void operator() ( size_type i, size_type j) {
		size_type siz = row_ptr(d.getSize());
		// if (bsp_pid() == 2)
		//	std::cout<<siz<<std::endl;
		// val(siz) = v;
		col_ind(siz) = j;
		size_type ind = i - d.getMinLocalGid();
		if (row_ptr(ind) == -1)
			row_ptr(ind) = siz;
		row_ptr(d.getSize()) ++;
	}

	Distribution const& distribution() const { return d ; }

	value_type getVal (size_type i) const {
		return 1;
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

        size_type num_rows() const {
		return nbrows_;	
	}

	size_type num_cols() const {
		return nbcols_;
	}

	};
}}

namespace glas2 {
template <typename D, typename T, typename F>
struct concept< bsp::BSP_COORD_CSR_matrix<D,T,F> >
: bsp::BSPCoordCsrMatrix
{};

}

