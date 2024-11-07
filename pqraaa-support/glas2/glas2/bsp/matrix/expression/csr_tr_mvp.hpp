//Zifan Liu, zifan.liu@cs.kuleuven.be

#ifndef glas2_bsp_matrix_expression_csr_tr_mvp_hpp

#define glas2_bsp_matrix_expression_csr_tr_mvp_hpp

#include <glas2/vector/concept/dense_vector.hpp>
//#include <glas2/vector/algorithm/inner_prod.hpp>
#include <cassert>
#include <iostream>

//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp.hpp"
//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp-templates.hpp"
#include <mcbsp.hpp>
#include <mcbsp-templates.hpp>


namespace glas2 { namespace bsp {

  template <typename M, typename V>
  class csr_tr_mvp {
    public:
      // value_type and size_type must be declared as a similar class to contiguous_vector class.
      typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;

    public:
      csr_tr_mvp( M const& m, V const& v )
      : m_( m )
      , v_( v )
      {
	value_type * z = new value_type[m_.distribution().getSize()];

	//First: compute the local MVP with sparse matrix in CSR format
	for (size_type i = 0; i < m_.distribution().getSize(); i++) {
		z[i] = 0;
		//if (bsp_pid()==0)
		//std::cout<<"enter here!!!!"<<std::endl;
		bsp_sync();
		for (size_type j = m_.getRow_ptr(i); j < m_.getRow_ptr(i+1); j++) {
			z[i] +=  m_.getVal(j) * v_(m_.getCol_ind(j)) ; 
		//	if (bsp_pid() == 0)
		//	std::cout<<"i = "<<i<<" j="<<j<<" row="<<m_.getRow_ptr(i)<<" col="<<m_.getCol_ind(j)<<" val="<<m_.getVal(j)<<std::endl;
		}
	}
	
	//Second step, locally computed output elements that should be stored at remote processes are sent out
	bsp_set_tagsize(sizeof(size_t));
	bsp_sync();
	for (size_type i = 0; i < bsp_nprocs(); ++i) {
		if (i != bsp_pid()) {
			// Tag
			size_t tag = static_cast<size_t>( bsp_pid() );
			// Send Data
			//bsp_send (i, &tag, z,  m_.distribution().getSize() * sizeof (value_type));
			bsp_send (i, &tag, z,  m_.distribution().getSize() );
		}
	}
	// Sync to ensure fan-out is done
	bsp_sync();

	// Collect remote contributions
	const void * msg_payload;
	const void * msg_tag;
	// First get the remote contribution
	vec = new value_type[m_.distribution().getNumGlobalElements()];
	while (bsp_hpmove(&msg_tag, &msg_payload) != SIZE_MAX) {
		const size_t remote_id = * static_cast<const size_t *>(msg_tag);
		const size_t len = static_cast<const size_t>( m_.distribution().getRemoteLength(remote_id) );
	        const size_t start = static_cast<const size_t>( m_.distribution().getRemoteStart(remote_id) );	
		const value_type * const dat = static_cast<const value_type*> (msg_payload);
		for (typename V::size_type i = 0; i < len; ++i)
			vec[start + i] = dat[i];
	}
	// Secondly, we store the local contribution
	const size_t start = m_.distribution().getMinLocalGid();
	const size_t len = m_.distribution().getSize();
	for (typename V::size_type i = 0; i < len; ++i)
		vec[start + i] = z[i];
	delete [] z;
	siz = m_.distribution().getNumGlobalElements();
      
      }



    public:
      size_type size() const { return siz ; }

      value_type operator() ( size_t i ) const { return vec[i] ; }

    public:
      M const& matrix() const { return m_ ; }
      V const& vector() const { return v_ ; }

    private:
      M m_ ;
      V v_ ;
      value_type * vec;//********

      size_type siz;
  } ;

} } // namespace glas2::bsp

namespace glas2 {

  template <typename M, typename V>
  struct concept< bsp::csr_tr_mvp<M,V> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
