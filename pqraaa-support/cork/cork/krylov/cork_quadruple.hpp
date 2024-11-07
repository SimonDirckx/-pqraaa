//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_krylov_cork_quadruple_hpp
#define cork_krylov_cork_quadruple_hpp

#include <cork/options/max_krylov_dimension.hpp>
#include <cork/options/value_of.hpp>
#include <cork/krylov/toar_triple.hpp>
#include <type_traits>
#include <cmath>

namespace CORK { namespace krylov {

  template <typename T, typename Options, typename PT=T>
  class cork_quadruple
  : public toar_triple<T, Options>
  {
    public:
      typedef toar_triple<T,Options> triple_type ;

      typedef typename triple_type::size_type    size_type ;
      typedef typename triple_type::value_type   value_type ;
      typedef typename triple_type::options_type options_type ;
      typedef typename triple_type::matrix_type  matrix_type ;
      typedef typename triple_type::sharedmatrix_type  sharedmatrix_type ;

      typedef PT pole_value_type ;

    public:
      cork_quadruple( size_type n, size_type k_max, size_type degree, options_type const& options )
      : triple_type( n, k_max, degree, options )
      //, continuation_combination( k_max, k_max )
      , K( k_max+1, k_max )
      {
        //glas2::fill( continuation_combination, 0.0 ) ;
      }

    public:
      template <typename Linearization>
      cork_quadruple( Linearization const& linearization, int n_wanted, options_type const& options )
      : cork_quadruple( linearization.size()
                      , (options::value_of<options::max_krylov_dimension>(options)==0? std::max(n_wanted+20,2*n_wanted) : options::value_of<options::max_krylov_dimension>(options))
                      , linearization.size_of_basis()
                      , options
                      )
      {}

      template <typename Linearization>
      cork_quadruple( Linearization const& linearization, options_type const& options )
      : cork_quadruple( linearization.size()
                      , options::value_of<options::max_krylov_dimension>(options)
                      , linearization.size_of_basis()
                      , options
                      )
      {
        assert( options::value_of<options::max_krylov_dimension>(options)!=0 ) ;
      }

    public:
      // Here the pole should have the same value_type as the quadruple.
      void next_step( size_type new_rank, value_type const pole, bool q_already_orto=true ) {
        static_cast<triple_type&>(*this).next_step( new_rank, q_already_orto ) ;
        auto column_K = K( glas2::all(), this->k()-1 ) ;
        column_K( glas2::range(0,this->k()+1) ) = this->H( glas2::range(0,this->k()+1), this->k()-1 ) * pole ;
        column_K( this->k()-1 ) += 1.0 ;
        fill( column_K(glas2::range_from_end(this->k()+1,0)), 0.0 ) ;
      } // next_step()

      // Is only possible for adding two complex conjugate solutions, i.e., real part and imaginary part in U.
      void next_step_double( size_type new_rank, pole_value_type const pole ) {
        auto column_K = K( glas2::all(), this->k() ) ;
        auto column_K1 = K( glas2::all(), this->k()+1 ) ;

        fill( column_K, 0.0 ) ;
        fill( column_K1, 0.0 ) ;

        static_cast<triple_type&>(*this).next_step( new_rank ) ;
        column_K( glas2::range(0,this->k()+1) ) = this->H( glas2::range(0,this->k()+1), this->k()-1 ) * real(pole) ;
        column_K( this->k()-1 ) += 1.0 ;
        column_K1( glas2::range(0,this->k()+1) ) = this->H( glas2::range(0,this->k()+1), this->k()-1 ) * imag(pole) ;

        static_cast<triple_type&>(*this).next_step( new_rank ) ;
        column_K( glas2::range(0,this->k()+1) ) -= this->H( glas2::range(0,this->k()+1), this->k()-1 ) * imag(pole) ;
        column_K1( glas2::range(0,this->k()+1) ) += this->H( glas2::range(0,this->k()+1), this->k()-1 ) * real(pole) ;
      } // next_step_double()

    public:
      //matrix_type continuation_combination ;
      sharedmatrix_type       K ; // Shared needed for invariant pair.
  } ; // class cork_quadruple
   
} } // namespace CORK::krylov


#endif
