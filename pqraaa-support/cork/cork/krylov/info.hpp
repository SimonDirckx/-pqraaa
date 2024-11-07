//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_info_hpp
#define cork_krylov_info_hpp

#include <cork/linearization/info.hpp>

namespace CORK { namespace krylov {

  struct info
  : public linearization::info
  {
    inline info()
    : krylov_process("NO INFORMATION")
    {}

    std::string                krylov_process ;
  } ;

  inline std::ostream& operator<<( std::ostream& os, info const& inf ) {
    os << static_cast<linearization::info const&>(inf) << "\n" ;
    os << "  type of Krylov method: " << inf.krylov_process << "\n" ;
    return os ;
  }
   
} } // namespace CORK::krylov


#endif
