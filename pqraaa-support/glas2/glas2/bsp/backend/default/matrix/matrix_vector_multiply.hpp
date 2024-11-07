// Zifan Liu, zifan.liu@cs.kuleuven.be

#ifndef glas2_bsp_backend_default_matrix_matrix_multiply_hpp
#define glas2_bsp_backend_default_matrix_matrix_multiply_hpp

#include <glas2/bsp/backend/default/default_backend.hpp>
#include <glas2/bsp/matrix/concept/bsp_row_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_column_matrix.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp.hpp"
//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp-templates.hpp"
#include <mcbsp.hpp>
#include <mcbsptemplates.hpp>
#include <glas2/matrix/algorithm/multiply.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>
//#include <glas2/bsp/matrix/expression/matrix_vector_multiply.hpp>
//#include <glas2/bsp/matrix/expression/matrix_multiply.hpp>

namespace glas2 { namespace bsp {
  

  template<typename M, typename V>
  typename std::enable_if< is<bsp::BSPRowMatrix, M>::value && is<bsp::BSPVector, V>::value 
  , bsp::distributed_vector< typename glas2::contiguous_vector<typename V::value_type,typename V::size_type>, typename M::distribution_type >
  >::type multiply(default_backend, M const& m, V const& v) {
		//First step, local MVP, since the portion of vector needed is local. No communication is needed here
		
	        typedef glas2::contiguous_vector<typename V::value_type,typename V::size_type> internal_type;
		typedef decltype( typename M::value_type() * typename V::value_type() ) value_type ;
		typedef decltype( typename M::size_type() + typename V::size_type() )   size_type ;
		//auto z = typename matrix_vector_multiply<typename M::local_type,V> multiply (m.local(), v);
		auto z = glas2::multiply (m.local(), v.local());

		//Second step, locally computed out- put elements that should be stored at remote processes are sent out
		bsp_set_tagsize(sizeof(size_t));
		size_t tag[2];

		for (typename V::size_type i = 0; i < bsp_nprocs(); ++i) {
			if (i != bsp_pid()) {
				// Tag
				tag[0] = static_cast<size_t>( m.distribution().getRemoteStart(i)  );
			        tag[1] = static_cast<size_t>( m.distribution().getRemoteLength(i)   );
				// Send Data
				bsp_hpsend (i, &tag, z +  tag[0], tag[1] * sizeof (value_type));
			}
		}
		// Sync to ensure fan-out is done
		bsp_sync();

		// Collect remote contributions
		const void * msg_payload;
		const void * msg_tag;
		while (bsp_hpmove(&msg_tag, &msg_payload) != SIZE_MAX) {
			const size_t start = static_cast<const size_t *>( msg_tag )[0];
			const size_t len = static_cast<const size_t*>( msg_tag )[1];
			const size_t * const dat = static_cast<const size_t*> (msg_payload);
			for (typename V::size_type i = 0; i < len; ++i)
				z[start + i] += dat[i]; 

		}
		const size_t start = v.distribution().getMinLocalGid(); 
		const size_t len = v.distribution().getSize();
		value_type * vec = new value_type[len];
		for (typename V::size_type i = 0; i < len; ++i)
			vec[i] = z[start + i];
		return bsp::distributed_vector< typename glas2::contiguous_vector<typename V::value_type,typename V::size_type>
			                                  , typename M::distribution_type
				>(vec, v.distribution());
	 }
	 

} } // namespace glas2::bsp

#endif
