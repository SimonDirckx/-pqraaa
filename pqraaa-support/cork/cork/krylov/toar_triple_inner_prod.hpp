//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_toar_triple_inner_prod_hpp
#define cork_krylov_toar_triple_inner_prod_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
//#include <boost/numeric/bindings/lapack/computational/ormqr.hpp>
//#include <boost/numeric/bindings/lapack/computational/unmqr.hpp>
//#include <boost/numeric/bindings/lapack/computational/geqrf.hpp>
//#include <boost/numeric/bindings/herm.hpp>
//#include <boost/numeric/bindings/trans.hpp>
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
  class toar_triple_inner_prod {
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
      toar_triple_inner_prod( size_type n, size_type k_max, size_type degree, options_type const& options )
      : Q( n, k_max+degree )
      , U( (k_max+degree)*degree, k_max+1 )
      , H( k_max+1, k_max )
      , rank( 0 )
      , k_max_( k_max )
      , k_( 0 )
      , degree_( degree )
      , options_( options )
      {
        assert( degree_ >= 2 ) ;
        glas2::fill( U, 0.0 ) ;
        glas2::fill( H, 0.0 ) ;
        //glas2::fill( Q, 0.0 ) ;
      }

    public:
      template <typename Linearization>
      toar_triple( Linearization const& linearization, int n_wanted, options_type const& options )
      : toar_triple( linearization.size()
                   , (options.maximum_krylov_dimension==0? 2*n_wanted : options.maximum_krylov_dimension)
                   , linearization.size_of_basis()
                   , options
                   )
      {}

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

      size_type k() const { return k_ ; }
      size_type degree() const { return degree_ ; }
      size_type k_max() const { return k_max_ ; }

    public:
      // step increases by 1.
      bool initial_vector( size_type new_rank, bool q_already_orto=true ) {
        rank = 0 ;
        k_ = 0 ;
        if (!q_already_orto) {
          for (size_type i=rank; i<new_rank; ++i) {
            add_column_of_Q( 0 ) ;
          }
        } else {
          rank = new_rank ;
        }
        auto nU = norm_2( U( glas2::all(), 0 ) ) ;
        if (nU!=0) {
          U( glas2::all(), 0 ) /= nU ;
        } else {
          return false ;
        }
        return true ;
      } // initial_vector()

      void next_step( size_type new_rank, bool q_already_orto=true ) {
        if (!q_already_orto) {
          for (size_type i=rank; i<new_rank; ++i) {
            add_column_of_Q( k_+1 ) ;
          }
        } else {
          rank = new_rank ;
        }
        gram_schmidt_u( k_+1 ) ;
        ++k_ ;
      } // next_step()

    private:
      void add_column_of_Q( size_type column ) {
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

      void gram_schmidt_u( size_type column ) {
        //auto norm_u = norm_2( U(glas2::all(),column) );

        glas2::vector< value_type > g( column ) ;
        glas2::range r_step(0,column) ;
        g = multiply( transpose(conj(U(glas2::all(),r_step))), U(glas2::all(),column) ) ;
        U(glas2::all(),column) -= multiply( U(glas2::all(),r_step), g ) ;
        H( r_step, column-1 ) = g ;
        g = multiply( transpose(conj(U(glas2::all(),r_step))), U(glas2::all(),column) ) ;
        U(glas2::all(),column) -= multiply( U(glas2::all(),r_step), g ) ;
        H( r_step, column-1 ) += g ;
        H( column, column-1 ) = norm_2( U(glas2::all(),column) ) ;

        fill( H(glas2::range_from_end(column+1,0), column-1), 0.0 ) ;

        // THIS IS STILL NOT SAFE:
        // WHEN ORTO FAILS: STOP OR CONTINUE WITH RANDOM VECTOR
        if (H(column,column-1)!=.0) U(glas2::all(),column) /= H(column,column-1) ;
        else
          std::cerr << "There is an upper Hessenberg matrix subdiagonal element zero" << std::endl ;
      } // gram_schmidt_u()

    public:
      // Multiply vectors on the right by P
      template <typename PP>
      void transform_vectors( PP const& P ) {
        assert( P.num_rows()==k()+1 ) ;
        assert( P.num_columns()<=k()+1 ) ;

        int block_size = k_max_+degree_ ;

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

        // Apply to Q
        /*        matrix_type Q_copy( Q.num_rows(), new_rank ) ;
        Q_copy = multiply( Q( glas2::all(), glas2::range(0,rank) ), orto_transform ) ;
        Q( glas2::all(), glas2::range(0,new_rank) ) = Q_copy ;
        */

        matrix_type UQ( rank, new_rank ) ;
        for (size_type i=0; i<Q.num_rows(); i+=rank) {
          size_type rQ = std::min(Q.num_rows()-i, rank) ;
          UQ( glas2::range(0,rQ), glas2::all() ) = multiply( Q( glas2::range(i, i+rQ), glas2::range(0,rank) ), orto_transform ) ;
          Q( glas2::range(i, i+rQ), glas2::range(0,new_rank) ) = UQ( glas2::range(0,rQ), glas2::all() ) ;
        }

        rank = new_rank ;
        k_ = new_k-1 ;
      } // transform()

    public:
      matrix_type Q ;
      matrix_type U ; // Why SHARED needed?
      matrix_type H ;
      size_type   rank ;

    private:
      size_type       k_max_ ;
      size_type       k_ ;
      size_type       degree_ ;
      options_type const& options_ ;
  } ; // class toar_triple
   
} } // namespace CORK::krylov


#endif
