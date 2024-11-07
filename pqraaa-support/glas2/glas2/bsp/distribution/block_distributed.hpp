//  (C) Copyright Karl Meerbergen & Albert-Jan Yzelman 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_vector_distribution_block_distributed_hpp
#define glas2_bsp_vector_distribution_block_distributed_hpp

#include <functional>
#include <type_traits>

namespace glas2 { namespace bsp {

  struct block_distributed {
    // Some info required about the block distribution

    bool operator==(block_distributed that ) const { return true; }

    template <typename D>
    typename std::enable_if< !std::is_same<D,block_distributed>::value, bool>::type operator==( D that ) const {
      return false ;
    }
  } ;

} } // namespace glas::bsp

#endif
