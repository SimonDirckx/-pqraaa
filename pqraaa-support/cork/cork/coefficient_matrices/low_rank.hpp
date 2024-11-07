//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_low_rank_hpp
#define cork_coefficient_matrices_low_rank_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  //
  // Square low_rank matrices.
  // Low rank matrices are represented by a CoefficientMatrices and a vector low_rank that indicated which columns of the matrices are nonzero.
  //
  template <typename CoefficientMatrices, typename LowRank>
  class low_rank
  {
    public:
      typedef typename std::decay<CoefficientMatrices>::type     matrices_type ;
      typedef typename std::decay<LowRank>::type                 low_rank_type ;
      typedef typename matrices_type::value_type                 value_type ;
      typedef typename matrices_type::size_type                  grade_type ;
      typedef typename matrices_type::size_type                  size_type ;

    public:
      low_rank( CoefficientMatrices coefficient_matrices, LowRank low_rank )
      : coefficient_matrices_( coefficient_matrices )
      : low_rank_( low_rank )
      {}

    public:
      size_type num_rows() const { return coefficient_matrices_.num_rows() ; }
      size_type num_columns() const { return coefficient_matrices_.num_columns() ; }

      grade_type num_matrices() const { return coefficient_matrices_.num_matrices() ; }

    public:
      auto operator() ( int i ) const {
        return coefficient_matrices_( i ) ;
      }

      template <typename M>
      void fill( int i, M m ) const {
        coefficient_matrices_.fill( i, m ) ;
      }

    public:
      low_rank_type const& low_rank() const{ return low_rank_ ; }

    public:
      template <typename X, typename W>
      void multiply_add( size_type i, X const& x, W& w ) const {
        coefficient_matrices_.multiply_add( i, x, w ) ;
      } // apply_scheduled()

      template <typename Coefs, typename Matrix>
      void accumulate( Coefs const& coefs, Matrix A ) const {
        coefficient_matrices_.accumulate( coefs, A ) ;
      } // accumulate()

    private:
      CoefficientMatrices               coefficient_matrices_ ;
      LowRank                           low_rank_ ;
  } ; // low_rank

  template <typename CoefficientMatrices, typename LowRank>
  low_rank<CoefficientMatrices, LowRank> make_low_rank_copy( CoefficientMatrices const& d, LowRank const& low_rank ) { return low_rank<CoefficientMatrices, LowRank>(d,low_rank); }

  template <typename CoefficientMatrices, typename LowRank>
  low_rank<CoefficientMatrices const&, LowRank const&> make_low_rank_lvalue_reference( CoefficientMatrices const& d, LowRank const& low_rank ) { return low_rank<CoefficientMatrices const&, LowRank const&>(d,low_rank); }

  template <typename CoefficientMatrices>
  struct is_low_rank
  : std::false_type
  {} ;

  template <typename CoefficientMatrices, typename LowRank>
  struct is_low_rank< low_rank<CoefficientMatrices, LowRank> >
  : std::true_type
  {} ;

} } // namespace CORK::coefficient_matrices

#endif
