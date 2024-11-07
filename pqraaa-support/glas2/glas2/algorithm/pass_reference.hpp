#ifndef glas2_algorithm_pass_reference_hpp
#define glas2_algorithm_pass_reference_hpp

namespace glas2 {

  template <typename E>
  E const& pass_reference( E const& e ) {
    return e ;
  }

} // namespace glas2

#endif
