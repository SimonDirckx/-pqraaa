//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_coefficient_matrices4cork_hpp
#define cork_coefficient_matrices_coefficient_matrices4cork_hpp

#include <cork/coefficient_matrices/coefficient_matrices4cork_range.hpp>
#include <cassert>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  // Square coefficient_matrices4cork matrices.
  template <typename ValueType, typename CoefficientMatrices>
  class coefficient_matrices4CORK
  {
    public:
      typedef ValueType                                value_type ;
      typedef typename CoefficientMatrices::size_type  size_type ;

    public:
      coefficient_matrices4CORK( CoefficientMatrices const& matrices )
      : accumulator_( 1, matrices.num_matrices() )
      , use_( matrices.num_matrices() )
      {}

    public:
      typedef glas2::shared_matrix< value_type > accumulator_type ; // Should be shared_matrix instead of matrix
      typename accumulator_type::base_type accumulator() { return accumulator_ ; }
      auto use() { return use_ ; }

    public:
      template <typename Range>
      coefficient_matrices4CORK_range< coefficient_matrices4CORK, Range > range( Range& r ) {
        return coefficient_matrices4CORK_range< coefficient_matrices4CORK, Range >( *this, r ) ;
      }

    public:
      accumulator_type const& accumulator() const { return accumulator_ ; }

    public:
      template <typename QQ>
      void initialize_schedule( QQ const& Q ) {
        if (Q.num_columns()!=accumulator_.num_rows()) accumulator_.resize( Q.num_columns(), accumulator_.num_columns() ) ;
        fill( accumulator_, 0.0 ) ;
        fill( use_, false ) ;
      } // initialize()

      template <typename QQ>
      void initialize_schedule( QQ const& Q, glas2::range const& r ) {
        initialize_schedule( Q ) ;
      } // initialize()

      template <typename Z>
      void schedule( size_type i, Z const& z ) {
        assert( z.size() == accumulator_.num_rows() ) ;
        if (norm_2(z)!=0.) {
          use_(i) = true ;
          accumulator_( glas2::all(), i ) += z ;
        }
      } // schedule()

    public:
      template <typename QQ, typename W>
      void apply_scheduled( CoefficientMatrices const& coef, QQ const& Q, W& w ) const {
        assert( w.size()==Q.num_rows() ) ;
        glas2::vector< value_type > temp( Q.num_rows() ) ;
        for (int i=0; i<accumulator_.num_columns(); ++i) {
          if (use_(i)) {
            temp = glas2::multiply( Q, this->accumulator_( glas2::all(),i) ) ;
            coef.multiply_add( i, temp, w ) ;
          }
        }
      } // apply_scheduled()

      template <typename W>
      void apply_scheduled( CoefficientMatrices const& coef, W& w ) const {
        assert( w.size()==accumulator_.num_rows() ) ;
        for (int i=0; i<accumulator_.num_columns(); ++i) {
          if (use_(i)) {
            coef.multiply_add( i, this->accumulator_( glas2::all(),i), w ) ;
          }
        }
      } // apply_scheduled()

    private:
      accumulator_type    accumulator_ ;
      glas2::shared_vector<bool> use_ ;
  } ; // coefficient_matrices4CORK


} } // namespace CORK::coefficient_matrices

#endif
