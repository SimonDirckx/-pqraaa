//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_barycentric_rational_real_hpp
#define cork_basis4cork_barycentric_rational_real_hpp

#include <cork/basis/barycentric_rational_real.hpp>
#include <cork/basis/explicit_matrices.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getrs.hpp>
#include <cassert>
#include <string>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename ValueType, typename Weights, typename Nodes, typename I>
  class basis4CORK< ValueType, basis::barycentric_rational_real<Weights,Nodes,I> >
  : public explicit_matrices< typename std::common_type< typename weights_type::value_type, typename nodes_type::value_type >::type >
  {
    public:
      typedef typename std::decay<Weights>::type weights_type ;
      typedef typename std::decay<Nodes>::type   nodes_type ;
      typedef ValueType                          value_type ;
      typedef I                                  size_type ;

    public:
      explicit basis4CORK( basis::barycentric_rational_real<Weights,Nodes,I> const& basis )
      : explicit_matrices< typename std::common_type< typename weights_type::value_type, typename nodes_type::value_type >::type >( basis.num_terms()-1 )
      , nodes_( basis.nodes() )
      {
        // Compute M and N
        auto M = this->M() ;
        auto N = this->N() ;
        fill(M, 0.0) ; M(0,0) = -1.0 ;
        fill(N, 0.0) ;

        int i_row = 1;
        int max_nodes = nodes_.size()-1;
        if (nodes_(max_nodes).imag()!=0) --max_nodes ;
        for (int i=0; i<max_nodes; ++i) {
          if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()==0.0) {
            M(0,i) = 1.0 ;

            M(i_row,i) = -weights_(i+1).real()*nodes_(i).real() ;
            N(i_row,i) = -weights_(i+1).real()*s_ ;

            M(i_row,i+1) = weights_(i).real()*nodes_(i+1).real() ;
            N(i_row,i+1) = weights_(i).real() ;
            ++i_row ;
          } else if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()!=0.0) {
            M(0,i) = 1.0 ;

            M(i_row,i) = -weights_(i+1).real()*nodes_(i).real() ;
            N(i_row,i) = -weights_(i+1).real() ;
            M(i_row+1,i) = -weights_(i+1).imag()*nodes_(i).real() ;
            N(i_row+1,i) = -weights_(i+1).imag() ;

            M(i_row,i+1) = weights_(i).real()*nodes_(i+1).real() ;
            N(i_row,i+1) = weights_(i).real() ;
            M(i_row+1,i+1) = weights_(i).real()*nodes_(i+1).imag() ;
            M(i_row,i+2) = -M(i_row+1,i+1) ;
            M(i_row+1,i+2) = M(i_row,i+1) ;
            N(i_row+1,i+2) = N(i_row,i+1) ;
            i_row += 2 ;
          } else if (nodes_(i).imag()!=0.0 && nodes_(i+2).imag()==0.0) {
            assert( 0 ) ;
          } else {
            M(0,i) = 2.0 ; M(0,i+1) = 0.0 ;

            auto w_n = nodes_(i)*weights_(i+2) ;
            lu_(i_row,i) = weights_(i+2).real()*s_ - real(w_n) ;
            lu_(i_row,i+1) = -weights_(i+2).imag()*s_ + imag(w_n) ;
            lu_(i_row+1,i) = -lu_(i_row,i+1) ;
            lu_(i_row+1,i+1) = lu_(i_row,i) ;

            w_n = nodes_(i+2)*weights_(i) ;
            lu_(i_row,i+2) = real(w_n) - weights_(i).real()*s_ ;
            lu_(i_row,i+3) = - imag(w_n) + weights_(i).imag()*s_ ;
            lu_(i_row+1,i+2) = -lu_(i_row,i+3) ;
            lu_(i_row+1,i+3) = lu_(i_row,i+2) ;

            ++i;
            i_row += 2 ;
          }
        }
        if (nodes_(nodes_.size()-1).imag()==0.0) {
          M(0,nodes_.size()-1) = 1.0 ;
        } else {
          M(0,nodes_.size()-2) = 2.0 ;
          M(0,nodes_.size()-1) = 0.0 ;
        }
      }

    public:
      // Basis4CORK
      int num_terms() const { return nodes_.size()+1 ; }
      int size() const { return num_terms() ; }

      typename std::add_lvalue_reference< typename std::add_const<Weights>::type >::type weights() const& {
        return weights_ ;
      }

      typename std::add_lvalue_reference< typename std::add_const<Nodes>::type >::type nodes() const& {
        return nodes_ ;
      }

    public:
      // Basis4CORK
      void shift( value_type s ) {
        s_ = s ;

        // Form original M-sN matrix
        fill( lu_(0,glas2::all()), 1.0 ) ;
        int i_row = 1;
        int max_nodes = nodes_.size()-1;
        if (nodes_(max_nodes).imag()!=0) --max_nodes ;
        for (int i=0; i<max_nodes; ++i) {
          if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()==0.0) {
            lu_(0,i) = 1.0 ;

            lu_(i_row,i) = weights_(i+1).real()*(s_-nodes_(i).real()) ;

            lu_(i_row,i+1) = weights_(i).real()*(nodes_(i+1).real()-s_) ;
            ++i_row ;
          } else if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()!=0.0) {
            lu_(0,i) = 1.0 ;

            lu_(i_row,i) = weights_(i+1).real()*(s_-nodes_(i).real()) ;
            lu_(i_row+1,i) = weights_(i+1).imag()*(s_-nodes_(i).real()) ;

            lu_(i_row,i+1) = weights_(i).real()*(nodes_(i+1).real()-s_) ;
            lu_(i_row+1,i+1) = weights_(i).real()*nodes_(i+1).imag() ;
            lu_(i_row,i+2) = -lu_(i_row+1,i+1) ;
            lu_(i_row+1,i+2) = lu_(i_row,i+1) ;
            i_row += 2 ;
          } else if (nodes_(i).imag()!=0.0 && nodes_(i+2).imag()==0.0) {
            assert( 0 ) ;
          } else {
            lu_(0,i) = 2.0 ; lu_(0,i+1) = 0.0 ;

            auto w_n = nodes_(i)*weights_(i+2) ;
            lu_(i_row,i) = weights_(i+2).real()*s_ - real(w_n) ;
            lu_(i_row,i+1) = -weights_(i+2).imag()*s_ + imag(w_n) ;
            lu_(i_row+1,i) = -lu_(i_row,i+1) ;
            lu_(i_row+1,i+1) = lu_(i_row,i) ;

            w_n = nodes_(i+2)*weights_(i) ;
            lu_(i_row,i+2) = real(w_n) - weights_(i).real()*s_ ;
            lu_(i_row,i+3) = - imag(w_n) + weights_(i).imag()*s_ ;
            lu_(i_row+1,i+2) = -lu_(i_row,i+3) ;
            lu_(i_row+1,i+3) = lu_(i_row,i+2) ;

            ++i;
            i_row += 2 ;
          }
        }
        if (nodes_(nodes_.size()-1).imag()==0.0) {
          lu_(0,nodes_.size()-1) = 1.0 ;
        } else {
          lu_(0,nodes_.size()-2) = 2.0 ;
          lu_(0,nodes_.size()-1) = 0.0 ;
        }
        std::cout << nice_print(lu_) << std::endl;
        assert( i_row==nodes_.size() ) ;

        int info = boost::numeric::bindings::lapack::getrf( lu_, pivots_ ) ;
        assert( info>=0 ) ;
        if (info!=0) throw std::string("STATE_SPACE: shift failed") ;
      }

      // Basis4CORK
      value_type shift() const { return s_ ; }

      // Basis4CORK
      template <typename UM, typename ZM>
      void multiply_N( UM const& U, ZM Z ) const {
        fill( Z(0,glas2::all()), 0.0 ) ;
        for (size_type i=1; i<weights_.size(); ++i) {
          if (weights_(i-1).imag()==0) {
            Z( i, glas2::all() ) = weights_(i-1) * U( i+1, glas2::all() ) - weights_(i) * U( i, glas2::all() ) ;
          } else {
            Z( i, glas2::all() ) = weights_(i-1) * U( i+1, glas2::all() ) - weights_(i) * U( i, glas2::all() ) ;
          }
        }
        int i_row = 1;
        for (int i=0; i<max_nodes; ++i) {
          if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()==0.0) {
            Z( i_row, glas2::all() ) = weights_(i) * U( i+2, glas2::all() ) - weights_(i+1) * U( i+1, glas2::all() ) ;
            ++i_row ;
          } else if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()!=0.0) {
            Z(i_row,glas2::all()) = -weights_(i+1).real()* U( i+1, glas2::all() ) + (weights_(i).real()*nodes_(i+1).real()) * U( i+2, glas2::all() ) ;
            Z(i_row+1,glas2::all()) = weights_(i).real()* U( i+3, glas2::all() )  - weights_(i+1).imag()* U( i+1, glas2::all() ) ;
            i_row += 2 ;
          } else if (nodes_(i).imag()!=0.0 && nodes_(i+2).imag()==0.0) {
            assert( 0 ) ;
          } else {
            Z(i_row,glas2::all()) = - weights_(i+2).real() * U( i+1, glas2::all() )  + weights_(i+2).imag() * U( i+2, glas2::all() )  + weights_(i).real() * U( i+3, glas2::all() )
                                - weights_(i).imag() * U( i+4, glas2::all() ) ;
            Z(i_row+1,glas2::all()) = - weights_(i+2).real() * U( i+2, glas2::all() ) - weights_(i+2).imag() * U( i+1, glas2::all() )
                                + weights_(i).imag() * U( i+3, glas2::all() ) + weights_(i).real() * U( i+4, glas2::all() ) ;
            ++i;
            i_row += 2 ;
          }
        }
      } // multiply_N()

      // Basis4CORK
      template <typename Z0, typename ZZ>
      void lower_solve_right_hand_side( Z0 const& z0, ZZ Z ) const {
        assert( Z.num_rows()==num_terms()-1 ) ;
        fill( Z, 0.0 ) ;
        Z( 0, glas2::all() ) = - z0 ;
      } // multiply_N()

      // Basis4CORK
      template <typename ZM>
      void solve( ZM Z ) const {
        assert( Z.num_rows()==num_terms()-1 ) ;
        glas2::matrix< typename ZM::value_type > Zc(Z.num_rows(), Z.num_columns() ) ;
        Zc = Z ;
        boost::numeric::bindings::lapack::getrs( lu_, pivots_, Zc ) ;
        Z = Zc ;
      } // solve()

    private:
      weights_type const& weights_ ;
      nodes_type const&   nodes_ ;
  } ; // barycentric_rational_real

} } // namespace CORK::basis4cork

#endif
