#ifndef cork_eigenvalue_selector_right_most_hpp
#define cork_eigenvalue_selector_right_most_hpp

#include <cork/vector.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <cmath>
#include <type_traits>

namespace CORK { namespace eigenvalue_selector {

  template <typename T>
  class right_most
  {
    public:
      typedef T                       value_type;
      typedef decltype(std::abs(T())) real_type ;
      typedef T                       shift_type ;

    private:
      static real_type constexpr xi = 0.5 ;

    public:
      right_most( int n_wanted, shift_type const& reference )
      : n_wanted_( n_wanted )
      , reference_( reference)
      , factor_( std::pow( (1.0+xi)/(1.0-xi), 2 ) )
      {}

    public:
      const int     n_wanted_max()     const { return n_wanted_;     }

    public:
      auto shifts( int n_shifts ) const {
        assert( n_shifts>=1 ) ;
        // Real shifts
        glas2::shared_vector<value_type> shift( n_shifts ) ;
        shift(0) = reference_ ;
        for (int i=1;i<shift.size();++i) shift(i) = factor_ * shift(i-1) ;
        return shift ;
      }

    public:
      // It is important that the first n_wanted indeces in the order vector
      // are the indices of the eigen vector that are wanted.
      template <typename Eigen, typename Order>
      int sort_eigenvalues( Eigen const& eigen, Order& order ) const {
        assert( order.size()==eigen.size() ) ;
        glas2::vector<real_type> distances( copy(-real(eigen)) ) ;
        order = glas2::range(0,eigen.size() ) ;
        glas2::sort(distances, order);
        return std::min<int>( eigen.size(), n_wanted_ ) ;
      } // sort_eigenvalues()

    private:
      int const  n_wanted_ ;
      shift_type reference_ ;
      real_type  factor_ ;
  } ; // class right_most

} } // namespace CORK::eigenvalue_selector

#endif
