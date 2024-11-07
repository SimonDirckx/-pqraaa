#ifndef glas2_type_copy_of_hpp
#define glas2_type_copy_of_hpp

#include <type_traits>

namespace glas2 {

  template <typename E>
  class copy_of {
    public:
      copy_of( E const& e )
      : e_( e )
      {}

      E const& expression() const { return e_ ; }

    private:
      E const& e_ ;
  } ;

} // namespace glas2

#endif
