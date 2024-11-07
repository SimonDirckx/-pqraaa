#ifndef glas2_type_pass_reference_hpp
#define glas2_type_pass_reference_hpp

namespace glas2 {

  template <typename E>
  struct pass_reference {
    static const bool pass_by_value = true ;

    typedef E type ;

    type const& operator() ( E const& e ) const {
      return e ;
    }
  } ;

} // namespace glas2

#endif
