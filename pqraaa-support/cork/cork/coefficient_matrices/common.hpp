//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_common_hpp
#define cork_coefficient_matrices_common_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <vector>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  // Square dense matrices.
  template <typename ValueType, typename SizeType>
  class common
  {
    public:
      typedef ValueType value_type ;
      typedef SizeType  grade_type ;

    public:
      common( grade_type num_matrices )
      : accumulator_( 1, num_matrices )
      {}

    public:
      template <typename QQ>
      void initialize_schedule( QQ const& Q, glas2::range const& r ) {
        assert( r.begin()>=0 && r.end()<=accumulator_.num_columns() ) ;
        range_ = r ;
        if (Q.num_columns()!=accumulator_.num_rows()) accumulator_.resize( Q.num_columns(), accumulator_.num_columns() ) ;
        fill( accumulator_, 0.0 ) ;
      } // initialize()

      template <typename Z>
      void schedule( degree_type i, Z const& z ) {
        assert( i>=range_.begin() && i<range_.end() ) ;
        assert( z.size() == acumulator_.num_rows() ) ;
        accumulator_( glas2::all(), i ) ;
      } // schedule()

    protected:
      glas2::range                range_ ;
      glas2::matrix< value_type > accumulator_ ;
  } ; // dense

} } // namespace CORK::coefficient_matrices

#endif
