//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_coefficient_matrices4cork_range_hpp
#define cork_coefficient_matrices_coefficient_matrices4cork_range_hpp

#include <cassert>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  // Square coefficient_matrices4cork matrices.
  template <typename CoefficientMatrices4CORK, typename Range>
  class coefficient_matrices4CORK_range
  {
    public:
      typedef CoefficientMatrices4CORK                           coefficient_matrix_4cork_type ;
      typedef typename coefficient_matrix_4cork_type::value_type value_type ;
      typedef typename coefficient_matrix_4cork_type::size_type  size_type ;
      typedef Range                                              range_type ;

    public:
      coefficient_matrices4CORK_range( coefficient_matrix_4cork_type& matrices, range_type const& range )
      : coefficient_matrix_4cork_( matrices )
      , range_( range )
      {}

    public:
      template <typename Z>
      void schedule( size_type i, Z const& z ) {
        assert( i>=0 && i<range_.size() ) ;
        coefficient_matrix_4cork_.schedule( range_(i), z ) ;
      } // schedule()

    private:
      coefficient_matrix_4cork_type&   coefficient_matrix_4cork_ ;
      range_type                       range_ ;
  } ; // coefficient_matrices4CORK


} } // namespace CORK::coefficient_matrices

#endif
