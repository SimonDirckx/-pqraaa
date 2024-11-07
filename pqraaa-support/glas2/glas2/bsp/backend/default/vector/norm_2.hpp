//Zifan Liu	 zifan.liu@gmail.com
#ifndef glas2_bsp_backend_default_vector_norm_hpp
#define glas2_bsp_backend_default_vector_norm_hpp

#include <glas2/bsp/backend/default/default_backend.hpp>
#include <glas2/vector/algorithm/norm_2.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>
#include <cmath>

//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp.hpp"
#include <mcbsp.hpp>
//#include "/Users/zil/glas/glas2/glas2/bsp/deps/bsp/mcbsp-templates.hpp"

namespace glas2 { namespace bsp {

  template <typename X>
  typename std::enable_if< is<BSPVector,X>::value
                         , decltype( std::abs(typename X::value_type() ))
                         >::type norm_2( default_backend, X const& x ) {

    typedef decltype( std::abs(typename X::value_type()) ) value_type; 
    //typedef double value_type;

    value_type sum = 0 ;
    
    
    //Calculate the local sum
    auto localSum = glas2::norm_2( x.local() ) ;
    

    //initialisation for BSP inner product
    //const size_t size = static_cast<size_t>( bsp_nprocs() * sizeof(value_type) );
    const size_t size = static_cast<size_t>( bsp_nprocs() );
    value_type *ip_buffer = new  value_type[bsp_nprocs()];
    bsp_push_reg( ip_buffer, size );
    bsp_sync();
    // BSP communication part
    for (unsigned int i = 0; i < bsp_nprocs() ; ++i)
	bsp_put(i, &localSum, ip_buffer, bsp_pid(), 1);
	    //bsp_put(i, &localSum, ip_buffer, bsp_pid() * sizeof(value_type), sizeof(value_type));
    bsp_sync();
    // Reduce to the final sum
    for (unsigned int i = 0; i < bsp_nprocs() ; ++i)
	    sum += std::pow(  ip_buffer[i] , 2 ) ; 
    
    bsp_pop_reg (ip_buffer);
    delete [] ip_buffer;

    // Return
    return std::sqrt( sum ) ;
  }
} } // namespace glas2

#endif
