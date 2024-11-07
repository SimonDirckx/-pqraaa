#ifndef glas2_concept_assert_hpp
#define glas2_concept_assert_hpp

#include <cassert>
#include <iostream>

namespace glas2 {
  void assert( bool condition, char* s) {
    if (!condition) {
      std::cerr << s << std::endl ;
      ::assert(condition)
    }
  }
}

#endif
