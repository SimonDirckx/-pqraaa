#ifndef cork_eigenvalue_selector_eigenvalue_selector_template_hpp
#define cork_eigenvalue_selector_eigenvalue_selector_template_hpp

#include <cassert>
#include <type_traits>
#include <glas2/vector.hpp>

namespace CORK { namespace eigenvalue_selector {

  template <typename T, typename Sort>
  class eigenvalue_selector_template
  {
    public:
      typedef T value_type;

    public:
      eigenvalue_selector_template( )
      {} // eigenvalue_selector_template()

    public:
     /**
      * Sort the eigenvalues and return how much should be used.
      *
      * @param eigen      The vector of eigenvalues that have to be sorted.
      *                   This vector should not be altered. 
      * @param order      A vector that contains the order of the sorted values
      *                   when the function returns.
      *
      * @return n_inside  Return how much eigenvalues should be used.
      */
      template <typename Eigen, typename Order>
      int sort_eigenvalues( Eigen const& eigen, Order& order ) const {
        static_assert( std::is_same< T, typename Eigen::value_type>::value, "CORK::
        for (int i=0; i<eigen.size(); i++) order(i) = i;
        return eigen.size();
      } // sort_eigenvalues()

    private:
      // class variables
      int dummy = 0;

  } ; // class eigenvalue_selector_template

} } // namespace CORK::eigenvalue_selector

#endif
