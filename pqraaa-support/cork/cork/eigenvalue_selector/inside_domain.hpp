#ifndef cork_eigenvalue_selector_inside_domain_hpp
#define cork_eigenvalue_selector_inside_domain_hpp

#include <cork/vector.hpp>
#include <cork/concept/domain.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <cmath>
#include <type_traits>

namespace CORK { namespace eigenvalue_selector {

  template <typename Domain>
#ifdef CORK_USE_CONCEPTS
    requires CORK::Domain< Domain >
#endif
  class inside_domain
  {
    public:
      typedef typename std::decay<Domain>::type::value_type value_type;
      typedef decltype(std::abs(value_type()))              real_type ;

    public:
      inside_domain( Domain const& domain, int n_wanted, int n_wanted_min, real_type const& tolerance=0. )
      : domain_( domain )
      , n_wanted_( n_wanted )
      , n_wanted_min_( n_wanted_min )
      , tolerance_( tolerance )
      {
        assert( n_wanted_>=n_wanted_min ) ;
      } // inside_domain()

      inside_domain( Domain const& domain, int n_wanted, real_type const& tolerance=0. )
      : domain_( domain )
      , n_wanted_( n_wanted )
      , n_wanted_min_( n_wanted_ )
      , tolerance_( tolerance )
      {} // inside_domain()

    public:
      Domain const& domain()       const { return domain_;       }
      const int     n_wanted_max()     const { return n_wanted_;     }

    public:
      auto shifts( int n_shifts ) const {
        return domain_.discretize_coarse( n_shifts ) ;
        //return std::move( domain_.discretize_coarse(n_shifts) ) ;
      }

    public:
      // It is important that the first n_wanted indeces in the order vector
      // are the indices of the eigen vector that are wanted.
      template <typename Eigen, typename Order>
      int sort_eigenvalues( Eigen const& eigen, Order& order ) const {
        static_assert( std::is_convertible< value_type, typename Eigen::value_type >::value, "CORK::eigenvalue_selector::inside_domain: Eigen does not have the value_type T" ) ;

        typedef decltype(std::abs(value_type())) real_type ;
        assert( order.size()==eigen.size() ) ;
        int n_wanted = 0;
        glas2::vector<real_type> distances(eigen.size());
        for (typename Eigen::size_type i=0; i<eigen.size(); i++) {
          distances(i) = domain_.distance(eigen(i));
        }
        order = glas2::range(0,eigen.size() ) ;
        glas2::sort(distances, order);
        int n_inside = 0 ;
        for (typename Eigen::size_type i=0; i<eigen.size(); i++) {
          if (distances(i) <= tolerance_) n_inside++;  
        }
        n_wanted = std::max(n_inside,n_wanted_min_) ;
        return n_wanted;
      } // sort_eigenvalues()

    private:
      Domain const& domain_ ;
      int const     n_wanted_ ;
      int const     n_wanted_min_ ;
      real_type     tolerance_ ;

  } ; // class inside_domain

} } // namespace CORK::eigenvalue_selector

#endif
