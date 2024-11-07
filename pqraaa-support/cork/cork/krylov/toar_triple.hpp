//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_toar_triple_hpp
#define cork_krylov_toar_triple_hpp

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
#include <cork/options/max_krylov_dimension.hpp>
#include <cork/options/backend.hpp>
#include <cork/options/q_drop_tol.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/value_of.hpp>
#include <cork/krylov/euclidean_inner_product.hpp>
#include <type_traits>
#include <limits>
#include <cmath>

namespace CORK { namespace krylov {

  //
  // Iteration vectors:
  //
  //
  // 
  template <typename T, typename Options, typename InnerProd=euclidean_inner_product>
  class toar_triple {
    public:
      typedef T                                                    value_type ;
      typedef glas2::shared_matrix< value_type >                   sharedmatrix_type ;
      typedef glas2::matrix< value_type >                          matrix_type ;
      typedef typename matrix_type::size_type                      size_type ;
      typedef InnerProd                                            inner_product_type ;
      
    private:
      typedef decltype(std::abs(value_type()))                     real_value_type ;

    public:
      typedef Options                                              options_type ;

    public:
      toar_triple( size_type n, size_type k_max, size_type degree, options_type const& options, inner_product_type const& inner_prod=inner_product_type() )
      : Q( n, k_max+degree )
      , U( (k_max+degree)*degree, k_max+1 )
      , H( k_max+1, k_max )
      , rank( 0 )
      , k_max_( k_max )
      , k_( 0 )
      , degree_( degree )
      , options_( options )
      , inner_prod_( inner_prod.generate( *this ) )
      {
        assert( degree_ >= 2 ) ;
        glas2::fill( U, 0.0 ) ;
        glas2::fill( H, 0.0 ) ;
        //glas2::fill( Q, 0.0 ) ;
      }

    public:
      template <typename Linearization>
      toar_triple( Linearization const& linearization, int n_wanted, options_type const& options, inner_product_type const& inner_prod=inner_product_type() )
      : toar_triple( linearization.size()
                   , (CORK::options::value_of<options::max_krylov_dimension>(options)==0? std::max(3*n_wanted, n_wanted+20) : CORK::options::value_of<options::max_krylov_dimension>(options))
                   , linearization.size_of_basis()
                   , options
                   , inner_prod
                   )
      {}

      template <typename Linearization>
      toar_triple( Linearization const& linearization, options_type const& options, inner_product_type const& inner_prod=inner_product_type() )
      : toar_triple( linearization.size()
                   , CORK::options::value_of<options::max_krylov_dimension>(options)
                   , linearization.size_of_basis()
                   , options
                   , inner_prod
                   )
      {
        assert( CORK::options::value_of<options::max_krylov_dimension>(options)!=0 ) ;
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

      size_type k() const { return k_ ; }
      size_type degree() const { return degree_ ; }
      size_type k_max() const { return k_max_ ; }

    public:
      // step increases by 1.
      bool initial_vector( size_type new_rank, bool q_already_orto=true ) {
        assert( new_rank<=Q.num_columns() ) ;
        auto U0 = u_vector(0) ;
        rank = 0 ;
        k_ = 0 ;
        if (!q_already_orto) {
          int last_i = new_rank ;
          for (size_type i=0; i<new_rank; ++i) {
            int old_rank = rank ;
            add_column_of_Q( 0 ) ;
            if (rank>old_rank) inner_prod_.add_column_of_Q( Q, i, CORK::options::value_of<options::backend_key>(options_) ) ;
            else {
              // Rank did not change: remove column of Q and column of U
              for (int j=rank; j<last_i; ++j) {
                Q(glas2::all(), j) = Q(glas2::all(), j+1) ;
                U0(glas2::all(), j) = U0(glas2::all(), j+1) ;
              }
              fill( U0(glas2::all(),last_i), 0.0) ;
              --last_i ;
            }
          }
        } else {
          for (size_type i=rank; i<new_rank; ++i) {
            inner_prod_.add_column_of_Q( Q, i, CORK::options::value_of<options::backend_key>(options_) ) ;
          }
          rank = new_rank ;
        }
        auto nU = inner_prod_.norm( u_vector(0)(glas2::all(),glas2::range(0,rank)) ) ;
        //auto nU = norm_2( U( glas2::all(), 0 ) ) ;
        if (nU!=0 && nU!=std::numeric_limits<real_value_type>::infinity()) {
          u_vector(0)(glas2::all(),glas2::range(0,rank)) /= nU ;
          //auto u0r = u_vector(0)(glas2::all(),glas2::range(0,rank)) ;
          //divides_assign( options_.backend_key, u0r, nU ) ;
        } else {
          return false ;
        }
        return true ;
      } // initial_vector()

      void back_track() {
        --k_ ;
      } // back_track()

      void next_step( size_type new_rank, bool q_already_orto=true ) {
        assert( new_rank<=Q.num_columns() ) ;
        if (!q_already_orto) {
          for (size_type i=rank; i<new_rank; ++i) {
            //int old_rank = rank ;
            add_column_of_Q( k_+1 ) ;
            //if (rank>old_rank) inner_prod_.add_column_of_Q( Q, i, CORK::options::value_of<options::backend_key>(options_) ) ;
          }
        }/* else {
          for (size_type i=rank; i<new_rank; ++i) {
            inner_prod_.add_column_of_Q( Q, i, CORK::options::value_of<options::backend_key>(options_) ) ;
          }
          assert( rank==new_rank) ;
          rank = new_rank ;
        } WAS NOT CORRECT*/
        gram_schmidt_u( k_+1 ) ;
        //std::cout << "U(0) " << u_vector_k(0) << std::endl ;
        ++k_ ;
      } // next_step()

    public:
      // Use inner product
      void add_column_of_Q( size_type column ) {
        auto u_v = u_vector( column ) ;

        auto norm_q = norm_2( CORK::options::value_of<options::backend_key>(options_), Q(glas2::all(),rank) ) ;
        //std::cout << "add " << Q(glas2::all(),rank) << std::endl ;

        glas2::shared_vector< value_type > g( rank ) ;
        glas2::range r_rank(0,rank) ;
        g = multiply( transpose(conj(Q(glas2::all(),r_rank))), Q(glas2::all(),rank) ) ;
        Q(glas2::all(),rank) -= multiply( Q(glas2::all(),r_rank), g ) ;
        u_v( glas2::all(), r_rank ) += outer_prod( u_v( glas2::all(), r_rank.size() ), g ) ;

        auto norm = norm_2( CORK::options::value_of<options::backend_key>(options_), Q(glas2::all(),rank) ) ;
        auto norm_1 = norm_q ;
        for (int count=0; count<2 && norm<0.707*norm_1; ++count) {
          g = multiply( transpose(conj(Q(glas2::all(),r_rank))), Q(glas2::all(),rank) ) ;
          Q(glas2::all(),rank) -= multiply( Q(glas2::all(),r_rank), g ) ;
          u_v( glas2::all(), r_rank ) += outer_prod( u_v( glas2::all(), r_rank.size() ), g ) ;
          norm_1 = norm ;
          norm = norm_2( CORK::options::value_of<options::backend_key>(options_), Q(glas2::all(),rank) ) ;
        }

        if (norm>CORK::options::value_of<options::q_drop_tol<real_value_type>>(options_)*norm_q) {
          Q(glas2::all(),rank) /= norm ;
          u_v( glas2::all(), rank ) *= norm ;
          inner_prod_.add_column_of_Q( Q, rank, CORK::options::value_of<options::backend_key>(options_) ) ;
          ++rank ;
        } else {
          fill( u_v( glas2::all(), rank ), 0.0 ) ;
        }
        //std::cout << "u_v1 :" << u_v( glas2::all(), glas2::range(0,rank) ) <<std::endl ;
        //fill( u_v( glas2::all(), glas2::range_from_end(rank,0) ), 0.0 ) ;
      } // add_column_of_Q()

    private:
      void gram_schmidt_u( size_type column ) {
        //auto norm_u = norm_2( U(glas2::all(),column) );

        glas2::vector< value_type > g( column ) ;
        glas2::range r_step(0,column) ;
        glas2::range r_rank(0,rank) ;
        //g = multiply( transpose(conj(U(glas2::all(),r_step))), U(glas2::all(),column) ) ;
        fill(H( glas2::all(), column-1 ), 0.0) ;

        real_value_type norm_u0 = inner_prod_.norm( u_vector(column)(glas2::all(),r_rank) ) ;
        if (norm_u0==0.0) throw std::runtime_error( "CORK: a vector added to the Krylov space is zero" ) ;

        real_value_type norm_u = norm_u0;
        for (int count=0; count<3 && (count==0 || norm_u<=0.707*norm_u0); ++count) {
          fill(g,0.0) ;
          for (int j=0; j<degree_; ++j) {
            inner_prod_( j, conj(u_block(j)( r_rank, r_step )), u_vector(column)(j, r_rank ), g ) ;
            //g += multiply( transpose(conj(u_block(j)( r_rank, r_step ))), u_vector(column)(j, r_rank ) ) ;
          }
          //U(glas2::all(),column) -= multiply( U(glas2::all(),r_step), g ) ;
          auto Uc = U(glas2::all(),column) ;
          minus_assign( CORK::options::value_of<options::backend_key>(options_), Uc, multiply( U(glas2::all(),r_step), g ) ) ;
          //H( r_step, column-1 ) += g ;
          H( r_step, column-1 ) += g ;

          norm_u0 = norm_u ;
          norm_u = inner_prod_.norm( u_vector(column)(glas2::all(),r_rank) ) ;
        }
        //g = multiply( transpose(conj(U(glas2::all(),r_step))), U(glas2::all(),column) ) ;

        //real_value_type norm_u = inner_prod_.norm( u_vector(column)(glas2::all(),r_rank) ) ;

/*        fill(g,0.0) ;
        for (int j=0; j<degree_; ++j) {
          inner_prod_( j, conj(u_block(j)( r_rank, r_step )), u_vector(column)(j, r_rank ), g ) ;
        }
        U(glas2::all(),column) -= multiply( U(glas2::all(),r_step), g ) ;
        H( r_step, column-1 ) += g ;
        //H( column, column-1 ) = norm_2( U(glas2::all(),column) ) ;
*/

        H( column, column-1 ) = norm_u ;

        // THIS IS STILL NOT SAFE:
        // WHEN ORTO FAILS: STOP OR CONTINUE WITH RANDOM VECTOR
        if (H(column,column-1)!=.0) {
          //U(glas2::all(),column) /= H(column,column-1) ;
          auto Uc = U(glas2::all(),column) ;
          divides_assign( CORK::options::value_of<options::backend_key>(options_), Uc, H(column,column-1) ) ;
        } else {
          std::cerr << "There is an upper Hessenberg matrix subdiagonal element zero" << std::endl ;
        }
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
        assert( new_rank<=Q.num_columns() ) ;
        // Truncate further if possible
        if (CORK::options::value_of<options::debug_level>(options_)>3) std::cout << "SVD of truncation of U " << svd_val << "  Drop tol = " << CORK::options::value_of<options::q_drop_tol<real_value_type>>(options_)*svd_val(1) << std::endl ;
        for (size_type i=1; i<new_rank; ++i) {
          if (svd_val(i)<CORK::options::value_of<options::q_drop_tol<real_value_type>>(options_)*svd_val(1)) {
            new_rank = i ;
            break ;
          }
        }
        if (CORK::options::value_of<options::debug_level>(options_)>2) std::cout << "rank change after basis transformation from " << rank << " to " << new_rank << std::endl ;

        auto orto_transform = U_copy(glas2::all(), glas2::range(0,new_rank)) ;

        auto UP2 = UP( glas2::range(0,new_rank), glas2::all() ) ;
        for (int i=0; i<degree_; ++i) {
          UP2 = multiply( transpose(conj(orto_transform)), U( glas2::range(i*block_size, i*block_size+rank), glas2::range(0,new_k) ) );
          fill( U( glas2::range(i*block_size,(i+1)*block_size), glas2::all()), 0.0 ) ;
          U( glas2::range(i*block_size, i*block_size+new_rank), glas2::range(0,new_k) ) = UP2 ;
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

        inner_prod_.transform_vectors( P ) ;
        //std::cout << "U(0) " << u_vector_k(0) << std::endl ;
      } // transform_vectors()

    public:
      decltype (auto) Q_k() const { return Q( glas2::all(), glas2::range(0,rank) ) ; }
      decltype(auto) u_vector_k( size_type i ) {
        return u_vector(i)( glas2::all(), glas2::range(0,rank) ) ;
      } // u_vector()
      decltype(auto) u_block_k( size_type j ) const { return u_block(j)( glas2::range(0,rank), glas2::range(0,k_) ) ; }
      decltype(auto) u_block_k1( size_type j ) const { return u_block(j)( glas2::range(0,rank), glas2::range(0,k_+1) ) ; }

    public:
      sharedmatrix_type Q ; // Shared needed for invariant pair
      sharedmatrix_type U ; // Why SHARED needed?
      sharedmatrix_type H ;
      size_type   rank ;

    private:
      size_type           k_max_ ;
      size_type           k_ ;
      size_type           degree_ ;
      options_type const& options_ ;
      inner_product_type  inner_prod_ ;
  } ; // class toar_triple


} } // namespace CORK::krylov


#endif
