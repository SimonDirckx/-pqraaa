//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_special_inner_product_hpp
#define cork_krylov_special_inner_product_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace krylov {

  //
  // op computes the innerproducy:
  //  auto op( DenseVector v, DenseVector w )
  //     v^T w (conjugation should not be done in op)
  //
  template <typename T, typename Op>
  class special_inner_product {
    public:
      inline special_inner_product( int start, Op const& op )
      : start_( start )
      , op_( op )
      , LQtQL_( 0, 0 )
      {}

      inline special_inner_product( int start, Op const& op, int n_Q, int degree )
      : start_( start )
      , op_( op )
      , LQtQL_( n_Q, n_Q )
      , degree_( degree )
      {}

    public:
      template <typename Triple>
      special_inner_product generate( Triple const& triple ) const {
        return special_inner_product( start_, op_, triple.Q.num_columns(), triple.degree() ) ;
      }

    public:
      template <typename QQ, typename Backend>
      void add_column_of_Q( QQ const& Q, int r, Backend const& backend ) {
        for (int i=0; i<r; ++i) {
          LQtQL_( i, r ) = op_( conj(Q( glas2::all(), i )), Q( glas2::all(), r ) ) ;
          LQtQL_( r, i ) = glas2::conj( LQtQL_( i, r ) ) ;
        }
        LQtQL_( r, r ) = op_( conj( Q( glas2::all(), r ) ) , Q( glas2::all(), r ) ) ;
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
        glas2::matrix<T> temp( v.num_columns(), v.num_rows()-start_ ) ;
        temp = multiply( LQtQL_( glas2::range(0,v.num_columns()), glas2::range(0,v.num_columns()) ), transpose(v(glas2::range_from_end(start_,0), glas2::all())) ) ;
        g += glas2::real( trace( multiply( conj(v(glas2::range_from_end(start_,0), glas2::all())), temp ) ) ) ;
        return std::sqrt(g) ;
      }

    private:
      int              start_ ;
      Op               op_ ;
      glas2::shared_matrix<T> LQtQL_ ;
      int              degree_ ;
  } ; // class special_inner_product

  template <typename T, typename Op>
  decltype (auto) make_special_inner_product( int start, Op const& op ) {
    return special_inner_product<T, Op>( start, op ) ;
  }
   
} } // namespace CORK::krylov


#endif
