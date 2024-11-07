//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_low_rank_inner_product_hpp
#define cork_krylov_low_rank_inner_product_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace krylov {

  //
  // Special inner product for low rank matrices in the matrix polynomial
  //
  template <typename T, typename Range, typename EnableIf=void>
  class low_rank_inner_product {
  } ;

  //
  // Let
  //
  //   P(s) = A_0 + ... + phi_p(s) * A_d + phi_p+1(s) * L * R_p+1 + ... + phi_p+d(s) * L * R_p+d
  // where the action y = L * x is the same as y = x(range) where range is an argument of the constructor
  // The value of start corresponds to p above.
  //
  template <typename T, typename Range>
  class low_rank_inner_product< T, Range, typename std::enable_if< glas2::is< glas2::Vector, Range >::value >::type > {
    public:
      inline low_rank_inner_product( int start, Range const& range )
      : start_( start )
      , range_( range )
      , LQtQL_( 0, 0 )
      {}

      inline low_rank_inner_product( int start, Range const& range, int n_Q, int degree )
      : start_( start )
      , range_( range )
      , LQtQL_( n_Q, n_Q )
      , degree_( degree )
      {}

    public:
      template <typename Triple>
      low_rank_inner_product generate( Triple const& triple ) const {
        return low_rank_inner_product( start_, range_, triple.Q.num_columns(), triple.degree() ) ;
      }

    public:
      template <typename QQ, typename Backend>
      void add_column_of_Q( QQ const& Q, int r, Backend const& ) {
        LQtQL_( glas2::range(0,r+1), r ) = multiply( transpose( conj( Q(range_, glas2::range(0,r+1)) ) ), Q(range_, r) ) ;
        LQtQL_( r, glas2::range(0,r) ) = conj( LQtQL_( glas2::range(0,r), r ) ) ;
      }

      template <typename PP>
      void transform_vectors( PP const& P ) {
        glas2::matrix< T > temp( P.num_rows(), P.num_columns() ) ;
        temp = multiply( LQtQL_( glas2::range(0,P.num_rows()), glas2::range(0,P.num_rows()) ), P ) ;
        LQtQL_( glas2::range(0,P.num_columns()), glas2::range(0,P.num_columns()) ) = multiply( transpose(conj(P)), temp ) ;
      }

      template <typename W, typename V, typename G>
      void operator() ( int j, W const& w, V const& v, G g ) const {
        if (j<start_) {
          g += multiply( transpose(w), v ) ;
        } else {
          glas2::vector<T> temp( v.size() ) ;
          temp = multiply( LQtQL_( glas2::range(0,v.size()), glas2::range(0,v.size()) ), v) ;
          g += multiply( transpose(w), temp ) ;
        }
      }

      template <typename V>
      decltype(auto) norm( V const& v ) const {
        auto g = norm_fro_squared( v(glas2::range(0,start_), glas2::all() ) ) ;
        std::cout<<"LQtQL_ = " << LQtQL_( glas2::range(0,v.num_columns()), glas2::range(0,v.num_columns()) ) << std::endl ;
        std::cout<<"g = " << g << std::endl ;
        glas2::matrix<T> temp( v.num_columns(), v.num_rows()-start_ ) ;
        temp = multiply( LQtQL_( glas2::range(0,v.num_columns()), glas2::range(0,v.num_columns()) ), transpose(v(glas2::range_from_end(start_,0), glas2::all())) ) ;
        g += glas2::real( trace( multiply( conj(v(glas2::range_from_end(start_,0), glas2::all())), temp ) ) ) ;
        std::cout<<"g = " << g << std::endl ;
        return std::sqrt(g) ;
      }

    private:
      int              start_ ;
      Range            range_ ;
      glas2::matrix<T> LQtQL_ ;
      int              degree_ ;
  } ; // class low_rank_inner_product

  template <typename T, typename Range>
  decltype (auto) make_low_rank_inner_product( int start, Range const& range ) {
    return low_rank_inner_product<T, Range>( start, range ) ;
  }
   
} } // namespace CORK::krylov


#endif
