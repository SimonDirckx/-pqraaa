//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_detail_implicit_restart_hpp
#define cork_eigs_detail_implicit_restart_hpp

#include <cork/krylov/cork_quadruple.hpp>
#include <cork/krylov/toar_triple.hpp>
#include <cork/krylov/expand_to_full.hpp>
#include <cork/lapack/qr.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/filter_infinite_eigenvalue.hpp>
#include <cork/options/value_of.hpp>
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace CORK { namespace eigs { namespace detail {

 template <typename Triple, typename Options, typename QMatrix>
 void implicit_qr( Triple& triple, Options const& options, QMatrix& QM ) {
   typedef typename Triple::value_type            value_type ;
   typedef decltype( std::abs(value_type()) )     real_type ;

   glas2::matrix< value_type > V( triple.Q.num_rows()*triple.degree(), triple.k()+1 ) ;
   expand_to_full( triple, V ) ;

   if (options::value_of<options::filter_infinite_eigenvalue>(options)>0) {
     glas2::shared_vector< value_type > temp_all( triple.k()+1 ) ;
     glas2::shared_matrix< value_type > QL( triple.k()+1, triple.k()+1 ) ;

     if (options::value_of<options::debug_level>(options)>1) std::cout << "Filtering infinite eigenvalue with QR step" << std::endl ;

     auto Q = QM( glas2::range(0,triple.k()+1), glas2::range(0,triple.k()+1) ) ;
     Q = glas2::eye( Q.num_rows(), Q.num_columns() ) ;
     for (int irestart=0; irestart<options::value_of<options::filter_infinite_eigenvalue>(options); ++irestart) {
       glas2::range range_k(0,triple.k()-irestart) ;
       glas2::range range_k1(0,triple.k()+1-irestart) ;

       auto H_k = triple.H( range_k1, range_k ) ;
       auto Q_k = QL(range_k1, range_k1) ;
       auto temp = temp_all( range_k ) ;

       // Compute the SVD for the condition number of R
       glas2::matrix<value_type> HR( range_k1.size(), range_k.size() ) ; HR = H_k ;
       glas2::vector<real_type> svd( range_k.size() ) ;
#ifndef NDEBUG
       int info = 
#endif
       boost::numeric::bindings::lapack::gesvd( 'N', 'N', HR, svd, HR, HR ) ;
       assert( info==0 ) ;
       std::cout << "Condition number of H = " << svd(0) / svd(range_k.size()-1)<< std::endl ;

       // Compute QR factorization of H
#ifndef NDEBUG
       info =
#endif
       lapack::qr( H_k, Q_k ) ;
       assert( info==0 ) ;

       // Apply Q on the right of H.
       // Compute H = H * Q row by row; (since H is upper triangular, this is feasible in-place.)
       for (typename decltype(H_k)::size_type i=0; i<range_k.size(); ++i) {
         glas2::range r_i(i,range_k.size()) ;
         temp = multiply( transpose(Q_k(r_i, range_k) ), H_k(i, r_i ) ) ;
         H_k(i, glas2::all()) = temp ;
       }
       fill( H_k(range_k.size(), glas2::all()), 0.0 ) ;
       //fill( H_k(glas2::all(),range_k.size()-1), 0.0 ) ;

       for (int i=0; i<Q.num_rows(); ++i) {
         temp_all(range_k1) = multiply( transpose(Q_k), Q(i,range_k1) ) ;
         Q(i,range_k) = temp ;
       }
     }

     // Apply Q to the Krylov vectors
     triple.transform_vectors( Q(glas2::range(0,triple.k()+1), glas2::range(0,triple.k()+1-options::value_of<options::filter_infinite_eigenvalue>(options))) ) ;
   }
   auto V1 = V( glas2::all(), glas2::range(0,triple.k()+1) ) ;
   expand_to_full( triple, V1 ) ;
 } // implicit_qr()

 template <typename Quadruple, typename Options, typename QMatrix>
 void implicit_qz( Quadruple& quad, Options const& options, QMatrix& Q ) {
   typedef typename Quadruple::value_type value_type ;

   if (options::value_of<options::filter_infinite_eigenvalue>(options)>0) {
     assert(options::value_of<options::filter_infinite_eigenvalue>(options)==1) ;
     if (options::value_of<options::debug_level>(options)>1) std::cout << "Filtering infinite eigenvalue with QZ step" << std::endl ;

     glas2::range range_k(0,quad.k()) ;
     glas2::range range_k1(0,quad.k()+1) ;
     glas2::shared_vector< value_type > temp( range_k1.size() ) ;

     auto H_k = quad.H( range_k1, range_k ) ;
     auto K_k = quad.K( range_k1, range_k ) ;
     auto Q_k = Q(range_k1, range_k1) ;

     // Compute QR factorization of H
#ifndef NDEBUG
     int info =
#endif
     lapack::qr( H_k, Q_k ) ;
     assert( info==0 ) ;

     // Compute Q^* K column by column
     for (int i=0; i<K_k.num_columns(); ++i) {
       temp = multiply( transpose(conj(Q_k)), K_k(glas2::all(),i) ) ;
       K_k(glas2::all(),i) = temp ;
     }

     // Compute Z
     // Store householder vector in last row of K since we do not need this any longer.
     glas2::shared_vector<value_type> householder( quad.k() ) ;
     /*auto*/ householder = K_k( quad.k(), glas2::all() ) ;
     householder( quad.k()-1 ) += norm_2( householder ) * householder( quad.k()-1 ) / std::abs(householder( quad.k()-1 )) ;
     auto norm_householder = norm_2_squared( householder ) ;

     // Apply Z on the right of H.
     auto H_k1 = quad.H( range_k, range_k ) ;
     auto tempk = temp(range_k) ;
     tempk = multiply( H_k1, conj(householder) ) ;
     tempk *= 2.0 / norm_householder ;
     H_k1 -= outer_prod( tempk, householder ) ;
     fill( H_k(range_k.size(), glas2::all()), 0.0 ) ;

     // Apply Z on the right of K.
     auto K_k1 = quad.K( range_k, range_k ) ;
     tempk = multiply( K_k1, conj(householder) ) ;
     tempk *= 2.0 / norm_householder ;
     K_k1 -= outer_prod( tempk, householder ) ;
     fill( K_k(range_k.size(), glas2::all()), 0.0 ) ;

     // Apply Q to the Krylov vectors
     quad.transform_vectors( Q_k(glas2::all(), range_k) ) ;
   }
 } // implicit_qz()

} } } // namespace CORK::eigs::detail


#endif
