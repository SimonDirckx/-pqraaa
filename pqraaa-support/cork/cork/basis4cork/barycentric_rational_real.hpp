//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_barycentric_rational_real_hpp
#define cork_basis4cork_barycentric_rational_real_hpp

#include <cork/utility/is_infinite.hpp>
#include <cork/basis/barycentric_rational_real.hpp>
#include <cork/basis4cork/barycentric_rational_real_matrices.hpp>
#include <cork/basis4cork/explicit_matrices.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/basis4cork/evaluate.hpp>
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

  template <typename Weights, typename Nodes, typename I>
  class basis4CORK< basis::barycentric_rational_real<Weights,Nodes,I> >
  : public explicit_matrices< decltype( std::abs(typename basis::barycentric_rational_real<Weights,Nodes,I>::data_value_type()) ) >
  {
    public:
      typedef typename std::decay<Weights>::type weights_type ;
      typedef typename std::decay<Nodes>::type   nodes_type ;
      typedef I                                  size_type ;

    public:
      explicit basis4CORK( basis::barycentric_rational_real<Weights,Nodes,I> const& basis )
      : explicit_matrices< decltype( std::abs(typename basis::barycentric_rational_real<Weights,Nodes,I>::data_value_type()) ) >( basis.num_terms() )
      , weights_( basis.weights() )
      , nodes_( basis.nodes() )
      {
        barycentric_rational_real_matrices( weights_, nodes_, this->M(), this->N() ) ;
/*        // Compute M and N
        auto M = this->M() ;
        auto N = this->N() ;
        fill(M, 0.0) ; M(0,0) = -1.0 ;
        fill(N, 0.0) ;

        int i_row = 1;
        for (int i=0; i_row<M.num_rows(); ++i) {
          if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()==0.0) {
            if (is_infinite(nodes_(i))) {
              N(0,i+1) = -1.0 ;
              M(i_row,i+1) = weights_(i+1).real() ;

              M(i_row,i+2) = weights_(i).real()*nodes_(i+1).real() ;
              N(i_row,i+2) = weights_(i).real() ;
              assert(!is_infinite(nodes_(i+1))) ;
            } else {
              M(0,i+1) = 1.0 ;

              M(i_row,i+1) = -weights_(i+1).real()*nodes_(i).real() ;
              N(i_row,i+1) = -weights_(i+1).real() ;

              if (is_infinite(nodes_(i+1))) {
                M(i_row,i+2) = -weights_(i).real() ;
              } else {
                M(i_row,i+2) = weights_(i).real()*nodes_(i+1).real() ;
                N(i_row,i+2) = weights_(i).real() ;
              }
            }
            ++i_row ;
          } else if (nodes_(i).imag()==0.0 && nodes_(i+1).imag()!=0.0) {
            assert( !is_infinite(nodes_(i+1)) ) ;
            if (is_infinite(nodes_(i))) {
              N(0,i+1) = -1.0 ;

              M(i_row,i+1) = 2.0*weights_(i+1).real() ;
              M(i_row+1,i+1) = -2.0*weights_(i+1).imag() ;
            } else {
              M(0,i+1) = 1.0 ;

              M(i_row,i+1) = -2.0*weights_(i+1).real()*nodes_(i).real() ;
              N(i_row,i+1) = -2.0*weights_(i+1).real() ;
              M(i_row+1,i+1) = 2.0*weights_(i+1).imag()*nodes_(i).real() ;
              N(i_row+1,i+1) = 2.0*weights_(i+1).imag() ;
            }

            M(i_row,i+2) = weights_(i).real()*nodes_(i+1).real() ;
            N(i_row,i+2) = weights_(i).real() ;
            M(i_row+1,i+2) = -weights_(i).real()*nodes_(i+1).imag() ;
            M(i_row,i+3) = -M(i_row+1,i+2) ;
            M(i_row+1,i+3) = M(i_row,i+2) ;
            N(i_row+1,i+3) = N(i_row,i+2) ;

            i_row += 2 ;
          } else / *if (nodes_(i).imag()!=0.0)* / {
            assert( nodes_(i).imag()!=0.0 ) ;
            assert( !is_infinite(nodes_(i)) ) ;
            if (i_row>=M.num_rows()-1) {
              M(0,i+1) = 1.0 ; M(0,i+2) = 0.0 ;

              auto w_n = nodes_(i+1)*weights_(i) ;
              M(i_row,i+1) = imag(w_n) ;
              N(i_row,i+1) = weights_(i).imag() ;
              M(i_row,i+2) = real(w_n) ;
              N(i_row,i+2) = weights_(i).real() ;
              ++i_row ;
            } else {
              assert( nodes_(i+2).imag()!=0.0 ) ;

              M(0,i+1) = 1.0 ; M(0,i+2) = 0.0 ;

              auto w_n = nodes_(i)*weights_(i+2) ;
              M(i_row,i+1) = - real(w_n) ;
              N(i_row,i+1) = - weights_(i+2).real() ;
              M(i_row,i+2) = -imag(w_n) ;
              N(i_row,i+2) = -weights_(i+2).imag() ;
              M(i_row+1,i+1) = -M(i_row,i+2) ;
              N(i_row+1,i+1) = -N(i_row,i+2) ;
              M(i_row+1,i+2) = M(i_row,i+1) ;
              N(i_row+1,i+2) = N(i_row,i+1) ;

              w_n = nodes_(i+2)*weights_(i) ;
              M(i_row,i+3) = real(w_n) ;
              N(i_row,i+3) = weights_(i).real() ;
              M(i_row,i+4) = imag(w_n) ;
              N(i_row,i+4) = weights_(i).imag() ;
              M(i_row+1,i+3) = -M(i_row,i+4) ;
              N(i_row+1,i+3) = -N(i_row,i+4) ;
              M(i_row+1,i+4) = M(i_row,i+3) ;
              N(i_row+1,i+4) = N(i_row,i+3) ;

              ++i;
              i_row += 2 ;
            }
          }
        }
        if (nodes_(nodes_.size()-1).imag()==0.0) {
          M(0,nodes_.size()) = 1.0 ;
        } else {
          M(0,nodes_.size()-1) = 1.0 ;
          M(0,nodes_.size()) = 0.0 ;
        }
  //std::cout << "M = " << nice_print(M) << std::endl ;
  //std::cout << "N = " << nice_print(N) << std::endl ;
        */
      }

    public:
      // Basis4CORK
      int num_terms() const { return nodes_.size()+1 ; }

      typename std::add_lvalue_reference< typename std::add_const<Weights>::type >::type weights() const& {
        return weights_ ;
      }

      typename std::add_lvalue_reference< typename std::add_const<Nodes>::type >::type nodes() const& {
        return nodes_ ;
      }

    public:
      template <typename Coefs>
      void evaluate( Coefs coefs ) const {
        CORK::basis4CORK::evaluate( *this, coefs ) ;
      }

    private:
      weights_type const& weights_ ;
      nodes_type const&   nodes_ ;
  } ; // barycentric_rational_real

} } // namespace CORK::basis4cork

#endif
