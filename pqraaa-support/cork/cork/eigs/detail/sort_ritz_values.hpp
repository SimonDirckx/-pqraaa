//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2015.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_detail_sort_ritz_values_hpp
#define cork_eigs_detail_sort_ritz_values_hpp

#include <cork/eigs/info.hpp>
#include <cork/eigs/nep_residual_norm.hpp>
#include <cork/options/value_of.hpp>
#include <cork/options/debug_level.hpp>
#include <cork/options/restart_reduction_factor.hpp>
#include <cork/options/number_of_stagnation_restarts.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <iosfwd>
#include <cmath>

namespace CORK { namespace eigs { namespace detail {

  template <typename EigSelector, typename Options>
  class sort_ritz_values {
    public:
      template <typename Subdiagonal>
      sort_ritz_values( EigSelector const& eigenvalue_selector, Options const& options, info& information, int& stagnate, Subdiagonal const& sub_diag )
      : eigenvalue_selector_( eigenvalue_selector )
      , options_( options )
      , n_converged_unsorted_( 0 )
      , information_( information )
      , not_splittable_( sub_diag.size()+1 )
      , restart_reduction_factor_( std::min( std::max<int>(1,options::value_of<options::restart_reduction_factor>(options_)), 90) )
      , stagnate_( stagnate )
      , last_number_converged_and_wanted_(information.number_converged_and_wanted)
      {
#ifndef NDEBUG
        fill( not_splittable_, -1 ) ;
#endif
        not_splittable_(0) = 0 ;
        // Identify when eigenvalues cannot be separated (only for real matrices)
        for (decltype(sub_diag.size()) i=0; i<sub_diag.size(); ++i) {
          if (sub_diag(i)!=0.0) {
            not_splittable_(i+1) = not_splittable_(i) ;
            ++i ;
          }
          if (i<sub_diag.size()) not_splittable_(i+1) = not_splittable_(i)+1 ;
        }
      }

      template <typename RitzValues, typename Order>
      int wanted_first( RitzValues const& lambda_values, Order& order ) {
#ifndef NDEBUG
        auto k = lambda_values.size() ;
#endif

        assert( k==order.size() ) ;

        // Sort eigenvalues depending on the selection criterion, eigenvalue_selector
        int n_wanted = eigenvalue_selector_.sort_eigenvalues( lambda_values, order ) ;
        n_wanted = lambda_values.size() ;
        if (options::value_of<options::debug_level>(options_)>2) {
          std::cout << "All Ritz values sorted: " << lambda_values(order) << std::endl ;
        }

        return n_wanted ;
      } // wanted_first()

      template <typename RitzValues, typename Order>
      void step_one( RitzValues& lambda_values, Order& order, int eig_wanted ) {
        auto k = lambda_values.size() ;

        assert( k==order.size() ) ;

        // Sort eigenvalues depending on the selection criterion, eigenvalue_selector
        if (options::value_of<options::debug_level>(options_)>2) std::cout << "Ritz values: " << lambda_values << std::endl ;
        n_wanted_ = std::min( eig_wanted, eigenvalue_selector_.sort_eigenvalues( lambda_values, order ) ) ;
        if (n_wanted_>k-4) {
          n_wanted_ = std::max<int>(0,int(k)-4) ;
          information_.warnings.push_back("CORK: Krylov space dimension is too small compared to number of wanted eigenvalues") ;
          std::cout << "CORK: Krylov space dimension is too small compared to number of wanted eigenvalues\n" ;
          std::cout << "CORK: The number of wanted eigenvalues is reduced to " << n_wanted_ << "." << std::endl ;
        }

        n_keep_ = std::max<int>( n_wanted_+2, std::ceil((restart_reduction_factor_ * k)/100) ) ;
        if (n_keep_>=k) {
          std::cout << "Restart does not reduced the subspace. Reduce the option restart_reduction_factor." << std::endl ;
           information_.warnings.push_back("Restart does not reduced the subspace. Reduce the option restart_reduction_factor.") ;
        }
        if (not_splittable_(order(n_keep_-1))==not_splittable_(order(n_keep_))) ++n_keep_ ;

        if (options::value_of<options::debug_level>(options_)>2) {
          std::cout << "All Ritz values sorted: " << lambda_values(order) << std::endl ;
        }
        //if (std::norm(lambda_values(order(n_keep_-1))-std::conj(lambda_values(order(n_keep_))))==0.0) ++n_keep_ ;
      } // step_one()
 

      template <typename Quadruple, typename EigenVec, typename SchurVec, typename RitzValues, typename Order,  typename Resid, typename StopCriterion, typename T>
      bool step_two( Quadruple const& quadruple, EigenVec const& eigen_vec, SchurVec const& schur_vec, RitzValues const& lambda_values, Order& order, Resid const& resid,
                     StopCriterion const& stop_criterion, T const& norm_H ) {
        assert( lambda_values.size()==n_keep_ ) ;
        assert( lambda_values.size()==order.size() ) ;
        typedef typename RitzValues::value_type  value_type ;
        typedef decltype(std::abs(value_type())) real_type ;
        typedef std::complex<real_type>          complex_type ;

        glas2::vector<real_type> resid_nep( lambda_values.size() ) ;
        glas2::shared_vector<int> internal_order( lambda_values.size() ) ;
        int n_wanted = eigenvalue_selector_.sort_eigenvalues( lambda_values, internal_order ) ;
        //n_wanted = lambda_values.size() ; // WHY?
        assert( n_wanted <= n_keep_);

        // some preparation calculations for eigenvector if necessary
        glas2::matrix< complex_type>          combinations(0,0);
        glas2::shared_vector< complex_type >  eigvec(0);
        if ( stop_criterion.eigenvector_needed() ) {
          fill( resid_nep, std::numeric_limits<real_type>::infinity() ) ;

          glas2::matrix<value_type> schur_vec_times_H(schur_vec.num_rows(), eigen_vec.num_rows());
          glas2::matrix<complex_type> schur_vec_times_eigen_vec(schur_vec.num_rows(), eigen_vec.num_columns());
          combinations.resize( quadruple.rank, eigen_vec.num_columns() );

          schur_vec_times_H = multiply( schur_vec 
                                      , quadruple.H( glas2::range(0,schur_vec.num_columns()), glas2::range(0,eigen_vec.num_rows())) );
          schur_vec_times_eigen_vec = multiply( schur_vec_times_H, eigen_vec);
          combinations =      multiply( quadruple.U( glas2::range(0,quadruple.rank), glas2::range(0,schur_vec_times_eigen_vec.num_rows()) )
                                      , schur_vec_times_eigen_vec );
        }
        
        information_.number_converged = 0 ;
        information_.number_converged_and_wanted = 0;
        for (int i=0; i<n_keep_; ++i) {
          int index = internal_order(i);
          // calculating eigenvector if necessary
          if ( stop_criterion.eigenvector_needed() ){
            eigvec.resize(quadruple.Q.num_rows());
            eigvec =  multiply( quadruple.Q(glas2::all(), glas2::range(0,quadruple.rank))
                              , combinations(glas2::all(), index) );
            eigvec /= norm_2(eigvec);
          }
          if (stop_criterion.test( lambda_values(index), eigvec, resid(index), resid_nep(index), norm_H )) {
            order( information_.number_converged ) = internal_order(i) ;
            ++information_.number_converged ;
            if (i < n_wanted) { ++information_.number_converged_and_wanted; }
          } 
        }
        if (options::value_of<options::number_of_stagnation_restarts>(options_)<=0) information_.number_converged_and_wanted = information_.number_converged ;

        if (options::value_of<options::debug_level>(options_)>0) {
          std::cout << "Ritz values to be kept: " << lambda_values(internal_order) << std::endl ;
          std::cout << "Associated residual norms: " << resid(internal_order) << std::endl ;
          if (stop_criterion.eigenvector_needed()) std::cout << "Associated residual norms NEP: " << resid_nep(internal_order) << std::endl ;
        }


        bool converged = information_.number_converged_and_wanted >= n_wanted_ ;
        if (options::value_of<options::debug_level>(options_)>0) {
          std::cout << "Number of converged Ritz values before restart: " << information_.number_converged << "\n" ;
          std::cout << "Number of converged wanted Ritz values        : " << information_.number_converged_and_wanted << std::endl ;
        }

        if (!converged && information_.number_converged_and_wanted<=last_number_converged_and_wanted_ ) {
        //if (converged || (information_.number_converged_and_wanted==last_number_converged_ && information_.number_converged<=last_n_converged_unsorted) ) {
          //converged = false ;
          if (options::value_of<options::number_of_stagnation_restarts>(options_)>0 && stagnate_>=options::value_of<options::number_of_stagnation_restarts>(options_)) {
            if (options::value_of<options::number_of_stagnation_restarts>(options_)>0) {
              std::cout << "CORK: Possible stagnation detected: computations were stopped" << std::endl ;
              information_.warnings.push_back("Possible stagnation detected: computations were stopped") ;
            }
            information_.number_converged_and_wanted = information_.number_converged ;
            last_number_converged_and_wanted_ = information_.number_converged_and_wanted ;
            converged = true ;
          }
          ++ stagnate_ ;
        } else {
          last_number_converged_and_wanted_ = information_.number_converged_and_wanted ;
          stagnate_ = 0 ;
        }

        if (options::value_of<options::debug_level>(options_)>2) {
          if (information_.number_converged) std::cout << "Converged Ritz values: " << lambda_values(order(glas2::range(0,information_.number_converged))) << std::endl ;
          std::cout << "All Ritz values to be kept: " << lambda_values(internal_order) << std::endl ;
        }

        return converged ;
      } // step_two()


      int keep() const {return n_keep_ ;}

      void keep( int new_value ) {n_keep_ = new_value ;}

      int wanted() const {return n_wanted_ ;}

    private:
      EigSelector const&  eigenvalue_selector_ ;
      Options const&      options_ ;
      int                 n_wanted_ ;
      int                 n_keep_ ;
      int                 n_converged_unsorted_ ;
      info&               information_ ;
      glas2::shared_vector<int > not_splittable_ ;
      int                 restart_reduction_factor_ ;
      int&                stagnate_ ;
      int                 last_number_converged_and_wanted_ ;
  } ; // sort_ritz_values


  // lambda_values.size()==triple.k()

  // lambda_values.size()==keep


   template <typename EigSelector, typename RitzValues, typename RitzVectors>
   int sort_ritz_values_final( EigSelector const& eigenvalue_selector, RitzValues ritz_values, RitzVectors ritz_vectors ) {
     assert( ritz_vectors.num_columns() == ritz_values.size() ) ;

     auto k = ritz_values.size() ;
     glas2::vector<int> order( k ) ;

     // Sort eigenvalues depending on the selection criterion, eigenvalue_selector
     int n_wanted = eigenvalue_selector.sort_eigenvalues( ritz_values, order ) ;

     // Sort Ritz values accordingly.
     {
       glas2::vector< typename RitzValues::value_type > ritz_values_copy( glas2::copy( ritz_values ) ) ;
       ritz_values = ritz_values_copy( order ) ;
     }
     {
       glas2::matrix< typename RitzVectors::value_type > temp( ritz_vectors.num_rows(), k ) ;
       temp = ritz_vectors ;
       ritz_vectors = temp( glas2::all(), order ) ;
     }

     return n_wanted ;
   } // sort_ritz_values_final()

} } } // namespace CORK::eigs::detail


#endif
