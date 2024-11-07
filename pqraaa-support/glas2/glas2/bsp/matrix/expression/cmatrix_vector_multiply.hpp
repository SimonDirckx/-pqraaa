#ifndef glas2_bsp_matrix_expression_cmatrix_vector_multiply_hpp
#define glas2_bsp_matrix_expression_cmatrix_vector_multiply_hpp

#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <cassert>
#include <iostream>

//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp.hpp"

//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp-templates.hpp"
#include <mcbsp.hpp>
#include <mcbsp-templates.hpp>



namespace glas2 { namespace bsp {

  template <typename M, typename V>
  class cmatrix_vector_multiply {
    public:
      // value_type and size_type must be declared as a similar class to contiguous_vector class.
      typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;

    public:
      cmatrix_vector_multiply( M const& m, V const& v )
      : m_( m )
      , v_( v )
      {
//	std::cout<<"coucou"<<std::endl;
        //assert( v.size()==m.num_columns() ) ;
	//First: compute the local MVP
	auto mvp = glas2::multiply (m.local(), v); //z is of type matrix_vector_multiply

	value_type * z = new value_type[m_.num_rows()];
        for (int i = 0; i < m_.num_rows(); ++i) {
		z[i] = mvp(i);
	}	
	/*if (bsp_pid()==0) {
		std::cout<<"here!!!!!!! number of rows in this processor is "<<m_.num_rows()<<std::endl;
		for (int i = 0; i < m_.num_rows(); ++i)
			std::cout<<z[i]<<std::endl;
		std::cout<<"ahlala!!!!!!!"<<std::endl;
	}
	*/
	//Second step, locally computed output elements that should be stored at remote processes are sent out
	bsp_set_tagsize(sizeof(size_t));
	bsp_sync();
	//size_t tag[2];
	//size_t * tag = new size_t[2];
	for (typename V::size_type i = 0; i < bsp_nprocs(); ++i) {
		if (i != bsp_pid()) {
			// Tag
			size_t tag = static_cast<size_t>( bsp_pid() );
			// Send Data
			//if (i == 0)
			//	std::cout<<"tag = "<< tag << "z[0] = "<< z[0]<<std::endl;



			//bsp_send (i, &tag, z,  m_.num_rows() * sizeof (value_type));
			bsp_send (i, &tag, z,  m_.num_rows());
		//	bsp_send (i, &tag[0], z +  tag[0], tag[1] * sizeof (value_type));
		//	if (bsp_pid() == 0)
		//		std::cout<<"tag = "<< tag << " , !!z[0] = "<< z[0]<<std::endl;
		}
	}
	// Sync to ensure fan-out is done
	bsp_sync();

	// Collect remote contributions
	const void * msg_payload;
	const void * msg_tag;
	//size_t * msg_tag;
	// First get the remote contribution
	vec = new value_type[m_.distribution().getNumGlobalElements()];
	while (bsp_hpmove(&msg_tag, &msg_payload) != SIZE_MAX) {
		//const size_t start = static_cast<const size_t *>(msg_tag)[0];
		const size_t remote_id = * static_cast<const size_t *>(msg_tag);
		//const size_t len = static_cast<const size_t *>(msg_tag)[1] ;
		//const size_t start = static_cast<const size_t>( m_.distribution().getMinLocalGid() );
		const size_t len = static_cast<const size_t>( m_.distribution().getRemoteLength(remote_id) );
	        const size_t start = static_cast<const size_t>( m_.distribution().getRemoteStart(remote_id) );	
		//auto start = msg_tag[0];
		//auto len = msg_tag[1];
		const value_type * const dat = static_cast<const value_type*> (msg_payload);
		// Allocate the memory space for vec[]
		for (typename V::size_type i = 0; i < len; ++i)
			vec[start + i] = dat[i];

		//bsp_sync();
	/*	if (bsp_pid() == 1){
			std::cout<<"start="<<start<<" and len="<<len<<" for the remote id "<<remote_id<<std::endl;
			for (typename V::size_type i = 0; i < len; ++i)
				std::cout<<vec[start+i]<<std::endl;
		}
	*/	//bsp_sync(); Please be very careful with bsp_sync(), here we are deadlocked
	}
	//std::cout<<"iter = "<<iter<<std::endl;
	// Secondly, we store the local contribution
	const size_t start = m_.distribution().getMinLocalGid();
	const size_t len = m_.distribution().getSize();
	//if (bsp_pid() == 2)
	//	std::cout<<"pid = "<<bsp_pid()<<" start = "<<start<<" and len = "<<len<<" z[0] = "<<z[0]<<std::endl;
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
  struct concept< bsp::cmatrix_vector_multiply<M,V> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
