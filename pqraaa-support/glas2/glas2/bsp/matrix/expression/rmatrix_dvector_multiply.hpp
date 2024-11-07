#ifndef glas2_bsp_matrix_expression_rmatrix_dvector_multiply_hpp
#define glas2_bsp_matrix_expression_rmatrix_dvector_multiply_hpp

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
  class rmatrix_dvector_multiply {
    public:
      // value_type and size_type must be declared as a similar class to contiguous_vector class.
      typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
      typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;

    public:
      rmatrix_dvector_multiply( M const& m, V const& v )
      : m_( m )
      , v_( v )
      {
//	std::cout<<"coucou"<<std::endl;
        //assert( v.size()==m.num_columns() ) ;
	//First: compute the local MVP
	auto mvp = glas2::multiply (m.local(), v.local()); //z is of type matrix_vector_multiply

	value_type * z = new value_type[m_.num_rows()];
        for (int i = 0; i < m_.num_rows(); ++i) {
		z[i] = mvp(i);
	}	
/*	if (bsp_pid()==0) {
		std::cout<<"here!!!!!!!"<<std::endl;
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
			size_t tag[2];
			tag[0] = static_cast<size_t>( m_.distribution().getRemoteStart(i)  );
			tag[1] = static_cast<size_t>( m_.distribution().getRemoteLength(i) );
			// Send Data
			//bsp_send (i, &tag[0], z +  tag[0], tag[1] * sizeof (value_type));
			bsp_send (i, &tag[0], z +  tag[0], tag[1] );

			//bsp_send (i, &tag[0], z +  tag[0], tag[1] * sizeof (value_type));
			//if (i == 0)
			//	std::cout<<"tag[0]="<<tag[0]<<" and tag[1]="<<tag[1]<<std::endl;
		}
	}
	// Sync to ensure fan-out is done
	bsp_sync();

	// Collect remote contributions
	const void * msg_payload;
	const void * msg_tag;
	//size_t * msg_tag;
	while (bsp_hpmove(&msg_tag, &msg_payload) != SIZE_MAX) {
		//const size_t start = static_cast<const size_t *>(msg_tag)[0];
		const size_t start = * static_cast<const size_t *>(msg_tag);
		//const size_t len = static_cast<const size_t *>(msg_tag)[1] ;
		//const size_t start = static_cast<const size_t>( m_.distribution().getMinLocalGid() );
		const size_t len = static_cast<const size_t>( m_.distribution().getSize() ); 
		//auto start = msg_tag[0];
		//auto len = msg_tag[1];
		const value_type * const dat = static_cast<const value_type*> (msg_payload);
		for (typename V::size_type i = 0; i < len; ++i)
			z[start + i] += dat[i];
	//	if (bsp_pid() == 0)
	//		std::cout<<"!!start="<<start<<" and len="<<len<<std::endl;
	}
	const size_t start = v.distribution().getMinLocalGid();
	const size_t len = v.distribution().getSize();
	vec = new value_type[len];
	for (typename V::size_type i = 0; i < len; ++i)
		vec[i] = z[start + i];
	delete [] z;
	siz = len;
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
  struct concept< bsp::rmatrix_dvector_multiply<M,V> >
  : DenseVector
  {} ;

} // namespace glas2

#endif
