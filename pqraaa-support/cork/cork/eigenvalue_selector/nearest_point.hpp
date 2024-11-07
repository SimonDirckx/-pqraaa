#ifndef cork_eigenvalue_selector_nearest_point_hpp
#define cork_eigenvalue_selector_nearest_point_hpp

#include <cassert>
#include <type_traits>
#include <cork/vector.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace eigenvalue_selector {

  template <typename T>
  class nearest_point
  {
    public:
      typedef T                                value_type;
      typedef decltype(std::abs(value_type())) real_type ;

    public:
      nearest_point( value_type point, int n_wanted )
        : point_( point )
        , n_wanted_( n_wanted )
      {} // nearest_domain()

    public:
      // point is implemented such that it can be changed
      value_type&       point()        const { return point_;       }
      const int         n_wanted_max()     const { return n_wanted_;     }

    public:
      auto shifts( int n_shifts ) const {
        glas2::static_vector<T,1> shifts ; shifts(0) = point_ ;
        return shifts ;
      }

    public:
      // It is important that the first n_wanted indeces in the order vector
      // are the indeces of the eigen vector that are wanted.
      template <typename Eigen, typename Order>
      int sort_eigenvalues( Eigen const& eigen, Order& order ) const {
        static_assert( std::is_convertible< value_type, typename Eigen::value_type >::value, "CORK::eigenvalue_selector::nearest_point`: Eigen does not have the value_type T" ) ;

        assert( order.size()==eigen.size() ) ;
        glas2::vector<real_type> distances(eigen.size());
        for (typename Eigen::size_type i=0; i<eigen.size(); i++) {
          distances(i) = std::abs(point_-eigen(i));
        }
        order = glas2::range(0,eigen.size() ) ;
        glas2::sort(distances, order);
        return std::min<int>( eigen.size(), n_wanted_ ) ;
      } // sort_eigenvalues()

    private:
      value_type    point_ ;
      const int     n_wanted_ ;

  } ; // class nearest_domain

} } // namespace CORK::eigenvalue_selector

#endif
