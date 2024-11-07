//Zifan Liu, zifan.liu@cs.kuleuven.be
// parallel sparse matrix M1 * sequential dense matrix M2

#ifndef glas2_bsp_matrix_expression_coord_csr_mmp_hpp

#define glas2_bsp_matrix_expression_coord_csr_mmp_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>

#include <cassert>
#include <iostream>

#include <mcbsp.hpp>
#include <mcbsp-templates.hpp>


namespace glas2 { namespace bsp {

  template <typename M1, typename M2>
  class coord_csr_mmp {
    public:
      // value_type and size_type must be declared as a similar class to contiguous_vector class.
      typedef decltype( typename M1::value_type() * typename M2::value_type() ) value_type ;
      typedef decltype( typename M1::size_type() + typename M2::size_type() )   size_type ;

    public:
      coord_csr_mmp( M1 const& m1, M2 const& m2 )
      : m1_( m1 )
      , m2_( m2 )
      {
	nbrows_ = m1_.distribution().getNumGlobalElements();
	nbcols_ = m2_.num_columns();
	mat = new value_type[nbrows_ * nbcols_];

	const size_t start = static_cast<const size_t>( m1_.distribution().getMinLocalGid());
	//First: compute the local MMP with first sparse matrix in CSR format
	for (size_type i = 0; i < m1_.distribution().getSize(); i++) {
		for (size_type bk = 0; bk < nbcols_; bk ++)
			mat[i * nbcols_ + start + bk] = 0.0;
if (m1_.getRow_ptr(i) >= 0)
		for (size_type j = m1_.getRow_ptr(i) ; j < m1_.getRow_ptr(i+1); j++) {
			for  (size_type bk = 0; bk < nbcols_; bk ++) 
				mat[i * nbcols_ + start + bk] +=  m2_(m1_.getCol_ind(j), bk) ; 
		}
	}
	
	//Second step, locally computed output elements that should be stored at remote processes are sent out
	bsp_push_reg (mat, nbrows_ * nbcols_);
	bsp_sync();
	

	for (size_type i = 0; i < bsp_nprocs(); i++) {
		if (i != bsp_pid()) {
			bsp_put (i, mat + start * nbcols_, mat, start * nbcols_, m1_.distribution().getSize());
		}
	}
	// Sync to ensure fan-out is done
	bsp_sync();

	bsp_pop_reg (mat);

	bsp_sync();
      }



    public:
      size_type num_rows() const { return nbrows_ ; }
  
      size_type num_columns() const { return nbcols_ ; }

      value_type operator() ( size_t i, size_t j ) const { return mat[i * nbcols_ + j] ; }



    public:
      M1 const& sparse() const { return m1_ ; }
      M2 const& dense() const { return m2_ ; }

    private:
      M1 m1_ ;
      M2 m2_ ;
      value_type * mat;//********

      size_type nbrows_;
      size_type nbcols_;
  } ;

} } // namespace glas2::bsp

namespace glas2 {

  template <typename M1, typename M2>
  struct concept< bsp::coord_csr_mmp<M1,M2> >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
