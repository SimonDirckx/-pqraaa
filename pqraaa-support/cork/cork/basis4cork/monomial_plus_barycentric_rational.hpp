//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_monomial_plus_barycentric_rational_hpp
#define cork_basis4cork_monomial_plus_barycentric_rational_hpp

#include <cork/basis/monomial_plus_barycentric_rational.hpp>
#include <cork/basis/monomial.hpp>
#include <cork/basis4cork/explicit_matrices.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <string>

namespace CORK { namespace basis4CORK {

  template <typename BarycentricRational>
  class basis4CORK< basis::monomial_plus_barycentric_rational<BarycentricRational> >
  : public explicit_matrices< typename basis::monomial_plus_barycentric_rational<BarycentricRational>::barycentric_rational_type::minimal_value_type >
  {
    private:
      typedef explicit_matrices< typename basis::monomial_plus_barycentric_rational<BarycentricRational>::barycentric_rational_type::minimal_value_type > explicit_matrices_type ;

    public:
      typedef basis::monomial_plus_barycentric_rational<BarycentricRational> basis_type ;

      template< typename T>
      using value_type_for = typename basis_type::template value_type_for<T> ;

    private:
      typedef typename basis::monomial_plus_barycentric_rational<BarycentricRational>::barycentric_rational_type::minimal_value_type minimal_value_type ;

    public:
      // The size is the number of terms minus two, because the terms corresponding to the degree zero and maximum degree of monomial part are removed.
      explicit basis4CORK( basis_type const& basis )
      : explicit_matrices< minimal_value_type >( basis.num_terms()-2 )
      {
        // Compute M and N
        auto M = this->M() ;
        auto N = this->N() ;
        fill(M, 0.0) ;
        fill(N, 0.0) ;

        auto const& monomial = basis.monomial() ;
        auto const& barycentric = basis.barycentric() ;

        glas2::matrix< minimal_value_type > M_bary( barycentric.num_terms()-1, barycentric.num_terms() ) ;
        glas2::matrix< minimal_value_type > N_bary( barycentric.num_terms()-1, barycentric.num_terms() ) ;

        // Barycentric part
        basis4CORK< typename basis_type::barycentric_rational_type > bary_4cork( barycentric ) ;
        bary_4cork.fill_M( M_bary ) ;
        bary_4cork.fill_N( N_bary ) ;


        if (monomial.grade()>1) {
         // Monomial part
          fill( diagonal(M,0)(glas2::range(0,monomial.grade()-1)), 1.0 ) ;
          fill( diagonal(N,-1)(glas2::range(0,monomial.grade()-2)), 1.0 ) ;

          // First row
          M(0,0) = 1. ;
          N( 0,glas2::range_from_end(monomial.grade()-1,0) ) = M_bary( 0, glas2::range_from_end(1,0)) ;

          M(glas2::range_from_end(monomial.grade()-1,0), glas2::range_from_end(monomial.grade()-1,0) ) = M_bary( glas2::range_from_end(1,0), glas2::range_from_end(1,0) ) ;
          N(glas2::range_from_end(monomial.grade()-1,0), glas2::range_from_end(monomial.grade()-1,0) ) = N_bary( glas2::range_from_end(1,0), glas2::range_from_end(1,0) ) ;
        } else {
          M = M_bary( glas2::range_from_end(1,0), glas2::range_from_end(1,0) ) ;
          N = N_bary( glas2::range_from_end(1,0), glas2::range_from_end(1,0) ) ;
        }
      }

    public:
      template <typename T>
      using handle_type = explicit_matrices_handle< T,minimal_value_type > ;

      template <typename T>
      handle_type<T> handle() const {
        return handle_type<T>( this->M(), this->N(), []( T const& s ) { return s ; } ) ;
      }

    public:
      // Basis4CORK
      int num_terms() const { return this->size()+2 ; }
  } ; // monomial_plus_barycentric_rational

} } // namespace CORK::basis4cork

#endif
