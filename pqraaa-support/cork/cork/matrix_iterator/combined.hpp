//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_combined_hpp
#define cork_matrix_iterator_combined_hpp

#include <cork/coefficient_matrices/combined.hpp>
#include <cork/matrix_iterator/matrix_times_q.hpp>

namespace CORK { namespace matrix_iterator {

  template <typename Combined, typename CoefficientMatrices, typename Combinations, typename QQ>
  class matrix_times_q< Combined, QQ, combined< CoefficientMatrices, Combinations> >
  : public combined< matrix_times_q< typename std::decay<CoefficientMatrices>::type const&, QQ>, typename std::decay<Combinations>::type >
  {
    private:
      typedef matrix_times_q< typename std::decay<CoefficientMatrices>::type const&, QQ> matrix_times_q_type ;
      typedef combined< matrix_times_q_type, typename std::decay<Combinations>::type >   base_type ;

    public:
      matrix_times_q( Combined coefficient_matrices, QQ Q )
      : base_type( (matrix_times_q_type( coefficient_matrices.coefficient_matrices(), Q )), coefficient_matrices.combinations() )
      {}
  } ; // matrix_times_q


} } // namespace CORK::matrix_iterator

#endif
