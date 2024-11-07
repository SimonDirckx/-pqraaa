#ifndef cork_eigenvalue_selector_user_defined_hpp
#define cork_eigenvalue_selector_user_defined_hpp

#include <cassert>
#include <type_traits>
#include <cork/vector.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace eigenvalue_selector {

  template <typename T, typename Sort>
  class user_defined
  {
    public:
      typedef T value_type ;

    public:
      user_defined( int n_wanted, Sort const& sort )
      : n_wanted_( n_wanted )
      , sort_( sort )
      {} // nearest_domain()

    public:
      // point is implemented such that it can be changed
      const int         n_wanted_max()     const { return n_wanted_;     }

    public:
      // It is important that the first n_wanted indeces in the order vector
      // are the indeces of the eigen vector that are wanted.
      template <typename Eigen, typename Order>
      int sort_eigenvalues( Eigen const& eigen, Order& order ) const {
        static_assert( std::is_same< T, typename Eigen::value_type >::value, "CORK::eigenvalue_selector::user_defined: Eigen does not have the value_type T" ) ;
        return sort_( eigen, order ) ;
      } // sort_eigenvalues()

    private:
      const int     n_wanted_ ;
      Sort const&   sort_ ;

  } ; // class nearest_domain

} } // namespace CORK::eigenvalue_selector

#endif
