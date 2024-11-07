//Zifan Liu	 zifan.liu@gmail.com
#ifndef glas2_bsp_backend_default_vector_inner_prod_hpp
#define glas2_bsp_backend_default_vector_inner_prod_hpp

#include <glas2/bsp/backend/default/default_backend.hpp>
#include <glas2/vector/algorithm/inner_prod.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp.hpp"
#include <mcbsp.hpp>
//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp-templates.hpp"

namespace glas2 { namespace bsp {

  template <typename X, typename Y>
  typename std::enable_if< is<BSPVector,X>::value && is<BSPVector,Y>::value
                         , decltype( typename X::value_type() * typename Y::value_type() )
                         >::type inner_prod( default_backend, X const& x, Y const& y ) {
    assert( x.distribution()==y.distribution() ) ;
    typedef  decltype( typename X::value_type() * typename Y::value_type() ) value_type;
    
    value_type localSum = glas2::inner_prod( x.local(), y.local() ) ;
    
    //initialisation for BSP inner product
    //const size_t size = static_cast<size_t>( bsp_nprocs()* sizeof(value_type)  );
    const size_t size = static_cast<size_t>( bsp_nprocs() );
    value_type *ip_buffer = new  value_type[bsp_nprocs()];
    bsp_push_reg( ip_buffer, size );
    bsp_sync();
    // BSP communication part
    for (typename X::size_type i = 0; i < bsp_nprocs() ; ++i)
	bsp_put(i, &localSum, ip_buffer, bsp_pid(), 1);
	    //bsp_put(i, &localSum, ip_buffer, bsp_pid() * sizeof(value_type), sizeof(value_type));
    bsp_sync();

    // Reduce to the final sum
    for (typename X::size_type i = 1; i < bsp_nprocs() ; ++i)
	    ip_buffer[0] += ip_buffer[i]; 
    auto sum = ip_buffer[0];
    bsp_pop_reg (ip_buffer);
    delete [] ip_buffer;
    
    // Return
    return sum ;
  }
} } // namespace glas2

#endif
