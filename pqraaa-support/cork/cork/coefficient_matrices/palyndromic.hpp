//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_palyndromic_hpp
#define cork_coefficient_matrices_palyndromic_hpp

#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  template <typename CoefficientMatrices, typename NumMatrices=int>
  class palyndromic
  {
    public:
      typedef typename std::decay< CoefficientMatrices >::type   coefficient_matrices_type ;
      typedef typename coefficient_matrices_type::value_type     value_type ;
      typedef NumMatrices                                        grade_type ;
      typedef typename coefficient_matrices_type::size_type      size_type ;

    public:
      palyndromic( CoefficientMatrices coefficient, combinations )
      : coefficient_matrices_( coefficient )
      , num_matrices_( num_matrices )
      {
        assert( num_matrices_ >= 2*coefficient_matrices_.num_matrices()-1 ) ;
      }

    public:
      size_type size() const { return coefficient_matrices_.size() ; }

      grade_type num_matrices() const { return num_matrices_ ; }

    public:
      coefficient_matrices_type const& coefficient_matrices() const { return coefficient_matrices_ ; }

    private:
      NumMatrice num_matrices_ ;
  } ; // palyndromic


} } // namespace CORK::coefficient_matrices

#endif
