//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_compact_vectors_hpp
#define cork_krylov_compact_vectors_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <cork/krylov/options.hpp>
#include <type_traits>
#include <limits>
#include <cmath>

namespace CORK { namespace krylov {

  //
  // Iteration vectors:
  //
  //
  // 
  template <typename T>
  class compact_vectors {
    public:
      typedef T                                                    value_type ;
      typedef glas2::shared_matrix< value_type >                          sharedmatrix_type ;
      typedef glas2::matrix< value_type >                          matrix_type ;
      typedef typename matrix_type::size_type                      size_type ;
      
    private:
      typedef decltype(std::abs(value_type()))                     real_value_type ;

    public:
      typedef options< real_value_type >                           options_type ;

    public:
      compact_vectors( size_type n, size_type k_max, size_type degree, options_type const& options )
      : Q( n, k_max+degree )
      , U( (k_max+degree)*degree, k_max+1 )
      , rank( 0 )
      , k_max_( k_max )
      , degree_( degree )
      , options_( options )
      {
        assert( degree_ >= 2 ) ;
        glas2::fill( U, 0.0 ) ;
        //glas2::fill( Q, 0.0 ) ;
      }

    public:
      // This is U^T needed for Basis4CORK
      decltype(auto) u_vector( size_type i ) {
        assert( i>=0 && i<=k_max_ ) ;
        return reshape( U(glas2::all(), i), degree_, k_max_+degree_, glas2::row_major() ) ;
      } // u_vector()

      // Select horizontal block of U.
      decltype(auto) u_block( size_type j ) const {
        assert( j>=0 && j<degree_ ) ;
        return U( glas2::range(j*(k_max_+degree_), (j+1)*(k_max_+degree_)), glas2::all() ) ;
      } // u_block()

      size_type degree() const { return degree_ ; }
      size_type k_max() const { return k_max_ ; }

    public:
      // add_vector: (I\otimes Q) U(:,column)
      // Make Q orthonormal and update U(:,column)
      void add_vector( size_type column ) {
        auto u_v = u_vector( column ) ;

        auto norm_q = norm_2( Q(glas2::all(),rank) ) ;

        glas2::shared_vector< value_type > g( rank ) ;
        glas2::range r_rank(0,rank) ;
        g = multiply( transpose(conj(Q(glas2::all(),r_rank))), Q(glas2::all(),rank) ) ;
        Q(glas2::all(),rank) -= multiply( Q(glas2::all(),r_rank), g ) ;
        u_v( glas2::all(), r_rank ) += outer_prod( u_v( glas2::all(), r_rank.size() ), g ) ;

        g = multiply( transpose(conj(Q(glas2::all(),r_rank))), Q(glas2::all(),rank) ) ;
        Q(glas2::all(),rank) -= multiply( Q(glas2::all(),r_rank), g ) ;
        u_v( glas2::all(), r_rank ) += outer_prod( u_v( glas2::all(), r_rank.size() ), g ) ;

        auto norm = norm_2( Q(glas2::all(),rank) ) ;
        if (norm>options_.Q_drop_tol*norm_q) {
          Q(glas2::all(),rank) /= norm ;
          u_v( glas2::all(), rank ) *= norm ;
          ++rank ;
        }
        fill( u_v( glas2::all(), glas2::range_from_end(rank,0) ), 0.0 ) ;
      } // add_column_of_Q()

    public:
      // Multiply vectors on the right by P
      template <typename PP>
      void transform_vectors( PP const& P ) {
        assert( P.num_rows()<=k_max()+1 ) ;
        assert( P.num_columns()<=k_max()+1 ) ;

        int block_size = k_max()+degree() ;

        size_type new_k = P.num_columns() ;
        matrix_type UP( rank, new_k ) ;
        glas2::matrix<value_type> U_copy( rank, new_k*degree_ ) ;

        // Multiply the blocks of U by P
        for (int i=0; i<degree_; ++i) {
          UP = multiply( U( glas2::range(i*block_size, i*block_size+rank), glas2::range(0,k()+1) ), P );
          U( glas2::range(i*block_size, i*block_size+rank), glas2::range(0,new_k) ) = UP ;
          U_copy( glas2::all(), glas2::range(i*new_k, (i+1)*new_k) ) = UP ;
        }

        // Then, if possible remove rows of U and adapt Q and rank
        glas2::matrix<value_type> svd_U( 1, 1 ) ;
        glas2::vector<real_value_type> svd_val( std::min(rank, new_k*degree_) ) ;
        boost::numeric::bindings::lapack::gesvd( 'O', 'N', U_copy, svd_val, svd_U, svd_U ) ;
        //std::cout << "SVD: " << svd_val <<std::endl;
        size_type new_rank = std::min( new_k-1+degree_,svd_val.size() ) ; // This is smaller than new_k*degree_
        // Truncate further if possible
        if (options_.debug_level>3) std::cout << "SVD of truncation of U " << svd_val << std::endl ;
        for (size_type i=1; i<new_rank; ++i) {
          if (svd_val(i)<options_.Q_drop_tol*svd_val(1)) {
            new_rank = i ;
            break ;
          }
        }
        if (options_.debug_level>2) std::cout << "rank change after basis transformation from " << rank << " to " << new_rank << std::endl ;

        auto orto_transform = U_copy(glas2::all(), glas2::range(0,new_rank)) ;

        auto UP2 = UP( glas2::range(0,new_rank), glas2::all() ) ;
        for (int i=0; i<degree_; ++i) {
          UP2 = multiply( transpose(conj(orto_transform)), U( glas2::range(i*block_size, i*block_size+rank), glas2::range(0,new_k) ) );
          U( glas2::range(i*block_size, i*block_size+new_rank), glas2::range(0,new_k) ) = UP2 ;
          fill( U( glas2::range(i*block_size+new_rank,(i+1)*block_size), glas2::range(0,new_k)), 0.0 ) ;
        }

        matrix_type UQ( rank, new_rank ) ;
        for (size_type i=0; i<Q.num_rows(); i+=rank) {
          size_type rQ = std::min(Q.num_rows()-i, rank) ;
          UQ( glas2::range(0,rQ), glas2::all() ) = multiply( Q( glas2::range(i, i+rQ), glas2::range(0,rank) ), orto_transform ) ;
          Q( glas2::range(i, i+rQ), glas2::range(0,new_rank) ) = UQ( glas2::range(0,rQ), glas2::all() ) ;
        }

        rank = new_rank ;
      } // transform_vectors()

    public:
      decltype (auto) Q_k() const { return Q( glas2::all(), glas2::range(0,rank) ) ; }
      decltype(auto) u_vector_k( size_type i ) {
        return u_vector(i)( glas2::all(), glas2::range(0,rank) ) ;
      } // u_vector()
      decltype(auto) u_block_k( size_type j ) const { return u_block(j)( glas2::range(0,rank), glas2::range(0,k_) ) ; }
      decltype(auto) u_block_k1( size_type j ) const { return u_block(j)( glas2::range(0,rank), glas2::range(0,k_+1) ) ; }

    public:
      matrix_type Q ;
      matrix_type U ; // Why SHARED needed?
      size_type   rank ;

    private:
      size_type           k_max_ ;
      size_type           degree_ ;
      options_type const& options_ ;
  } ; // class compact_vectors


} } // namespace CORK::krylov


#endif
