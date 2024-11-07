//Zifan Liu, zifan.liu@cs.kuleuven.be

#ifndef glas2_bsp_matrix_expression_coord_csr_mvp_hpp

#define glas2_bsp_matrix_expression_coord_csr_mvp_hpp

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
  class coord_csr_mvp {
    public:
      // value_type and size_type must be declared as a similar class to contiguous_vector class.
      typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;

    public:
      coord_csr_mvp( M const& m, V const& v )
      : m_( m )
      , v_( v )
      {
	//value_type * z = new value_type[m_.distribution().getSize()];
	vec = new value_type[m_.distribution().getNumGlobalElements()];
	const size_t start = static_cast<const size_t>( m_.distribution().getMinLocalGid());
	//First: compute the local MVP with sparse matrix in CSR format
/*	for (unsigned int i = 0; i < bsp_nprocs(); i++){
		if (i == bsp_pid())
			std::cout<<"proc "<<i<<": size="<<m_.distribution().getSize()<<", start="<<start<<", global size="<<m_.distribution().getNumGlobalElements()<<std::endl;		
		bsp_sync();
	}
*/
	for (size_type i = 0; i < m_.distribution().getSize(); i++) {
		vec[i+start] = 0.0;
if (m_.getRow_ptr(i) >= 0)
		for (size_type j = m_.getRow_ptr(i) ; j < m_.getRow_ptr(i+1); j++) {
			//if ( m_.getCol_ind(j) >= 5563968) {
			//	std::cerr <<"fatal error vector size, proc="<<bsp_pid()<<" i="<<i<<" j="<<j<<" size="<<m_.getCol_ind(j)<<std::endl;
			//}
			
			vec[i+start] +=  v_(m_.getCol_ind(j)) ; 
		}
	}
	
	//Second step, locally computed output elements that should be stored at remote processes are sent out
	bsp_push_reg (vec, m_.distribution().getNumGlobalElements());
	bsp_sync();
	

	for (size_type i = 0; i < bsp_nprocs(); i++) {
		if (i != bsp_pid()) {
			bsp_put (i, vec+start, vec, start, m_.distribution().getSize());
		}
	}
	// Sync to ensure fan-out is done
	bsp_sync();

	bsp_pop_reg (vec);

	bsp_sync();
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
  struct concept< bsp::coord_csr_mvp<M,V> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
