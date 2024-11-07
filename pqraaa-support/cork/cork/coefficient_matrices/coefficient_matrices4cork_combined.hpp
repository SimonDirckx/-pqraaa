//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_coefficient_matrices4cork_combined_hpp
#define cork_coefficient_matrices_coefficient_matrices4cork_combined_hpp

#include <cork/coefficient_matrices/coefficient_matrices4cork_range.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  // Square coefficient_matrices4cork matrices.
  template <typename ValueType, typename CoefficientMatrices, typename Combinations>
  class coefficient_matrices4CORK< ValueType, combined<CoefficientMatrices,Combinations> >
  {
    public:
      typedef ValueType                                        value_type ;
      typedef typename std::decay< CoefficientMatrices >::type coefficient_matrices_type ;
      typedef typename coefficient_matrices_type::size_type    size_type ;

    private:
      typedef typename combined<CoefficientMatrices,Combinations>::combinations_type combinations_type ;

    public:
      coefficient_matrices4CORK( combined<CoefficientMatrices,Combinations> const& matrices )
      : accumulator_( 1, matrices.num_matrices() )
      , inner_accumulator_( 1, matrices.coefficient_matrices().num_matrices() )
      , combinations_( matrices.combinations() )
      {}

    public:
      typedef glas2::shared_matrix< value_type > accumulator_type ; // Should be shared_matrix instead of matrix
      typename accumulator_type::base_type accumulator() { return accumulator_ ; }

    public:
      template <typename Range>
      coefficient_matrices4CORK_range< coefficient_matrices4CORK, Range > range( Range& r ) {
        return coefficient_matrices4CORK_range< coefficient_matrices4CORK, Range >( *this, r ) ;
      }

    public:
      template <typename QQ>
      void initialize_schedule( QQ const& Q ) {
        if (Q.num_columns()!=accumulator_.num_rows()) {
          accumulator_.resize( Q.num_columns(), accumulator_.num_columns() ) ;
          inner_accumulator_.resize( Q.num_columns(), inner_accumulator_.num_columns() ) ;
        }
        fill( accumulator_, 0.0 ) ;
      } // initialize()

      template <typename QQ>
      void initialize_schedule( QQ const& Q, glas2::range const& ) {
        initialize_schedule( Q ) ;
      } // initialize()

      template <typename Z>
      void schedule( size_type i, Z const& z ) {
        assert( i>=0 && i<combinations_.num_columns() ) ;
        assert( z.size() == accumulator_.num_rows() ) ;
        accumulator_( glas2::all(), i ) += z ;
        //accumulator_( glas2::all(), i ) += multiply( combinations_, z ) ;
      } // schedule()

    public:
      template <typename QQ, typename W>
      void apply_scheduled( combined<CoefficientMatrices,Combinations> const& coef, QQ const& Q, W& w ) const {
        assert( w.size()==Q.num_rows() ) ;
        glas2::vector< value_type > temp( Q.num_rows() ) ;
        inner_accumulator_ = multiply( accumulator_, transpose(combinations_) ) ;
        for (int i=0; i<inner_accumulator_.num_columns(); ++i) {
          if (norm_2(this->inner_accumulator_( glas2::all(),i))!=0.) {
            temp = glas2::multiply( Q, this->inner_accumulator_( glas2::all(),i) ) ;
            coef.coefficient_matrices().multiply_add( i, temp, w ) ;
          }
        }
      } // apply_scheduled()

      template <typename W>
      void apply_scheduled( combined<CoefficientMatrices,Combinations> const& coef, W& w ) const {
        inner_accumulator_ = multiply( accumulator_, transpose(combinations_) ) ;
        for (int i=0; i<inner_accumulator_.num_columns(); ++i) {
          if (norm_2(this->inner_accumulator_( glas2::all(),i))!=0.) {
            coef.coefficient_matrices().multiply_add( i, this->inner_accumulator_( glas2::all(),i), w ) ;
          }
        }
      } // apply_scheduled()

    private:
      accumulator_type         accumulator_ ;
      mutable accumulator_type inner_accumulator_ ;
      combinations_type const& combinations_ ;
  } ; // coefficient_matrices4CORK

} } // namespace CORK::coefficient_matrices

#endif
