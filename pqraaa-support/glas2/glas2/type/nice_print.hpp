#ifndef glas2_type_nice_print_hpp
#define glas2_type_nice_print_hpp

#include <iosfwd>

namespace glas2 {

  
  template <typename E>
  struct nice_print_type {
    nice_print_type( E const& e )
    : e_( e )
    {}

    E const& e_ ;
  } ;

} // namespace glas2

#endif
