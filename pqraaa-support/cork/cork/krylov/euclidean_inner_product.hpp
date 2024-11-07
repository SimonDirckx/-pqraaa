//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_euclidean_inner_product_hpp
#define cork_krylov_euclidean_inner_product_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>

namespace CORK { namespace krylov {

  class euclidean_inner_product {
    public:
      inline euclidean_inner_product()
      {}

    public:
      template <typename Triple>
      euclidean_inner_product generate( Triple const& ) const { return *this ; }

    public:
      template <typename QQ, typename Backend>
      void add_column_of_Q( QQ const& Q, int r, Backend const& backend ) {
      }

      template <typename PP>
      void transform_vectors( PP const& P ) {
      }

      template <typename W, typename V, typename G>
      void operator() ( int, W const& w, V const& v, G g ) const {
        g += multiply( transpose( w ), v ) ;
      }

      template <typename V>
      decltype(auto) norm( V const& v ) const {
        return norm_fro( v ) ;
      }
  } ; // class euclidean_inner_product
   
} } // namespace CORK::krylov


#endif
