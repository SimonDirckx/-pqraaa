//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_selection_hpp
#define cork_coefficient_matrices_selection_hpp

#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  //
  // This class is used in combination with matrix_iterator< T, basis::merged >
  //
  template <typename CoefficientMatrices, typename Selection>
  class selection
  {
    public:
      typedef typename std::decay< CoefficientMatrices >::type                                                                         coefficient_matrices_type ;
      typedef typename std::decay< Selection >::type                                                                                selection_type ;
      typedef typename std::common_type< typename coefficient_matrices_type::value_type, typename selection_type::value_type>::type value_type ;
      typedef typename coefficient_matrices_type::grade_type                                                                           grade_type ;
      typedef typename coefficient_matrices_type::size_type                                                                            size_type ;

    public:
      selection( CoefficientMatrices coefficient, Selection selection )
      : coefficient_matrices_( coefficient )
      , selection_( selection )
      {}

    public:
      size_type num_rows() const { return coefficient_matrices_.num_rows() ; }
      size_type num_columns() const { return coefficient_matrices_.num_columns() ; }

      grade_type num_matrices() const { return selection_.size() ; }

    public:
      coefficient_matrices_type const& coefficient_matrices() const { return coefficient_matrices_ ; }
      selection_type const& my_selection() const { return selection_ ; }

    public:
      // Not efficient for a range of i's.
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        assert( i>=0 && i<num_matrices() ) ;
        assert( selection_(i)>=0) ;
        assert( selection_(i)<coefficient_matrices_.num_natrices()) ;
        coefficient_matrices_.multiply_add( selection_(i), x, w ) ;
      }

      template <typename Range, typename Lambda, typename W>
      void accumulate( Range const& range, Lambda const& lambda, W w ) const {
        coefficient_matrices_.accumulate( selection_(range), lambda, w ) ;
      } // accumulate()

    private:
      CoefficientMatrices coefficient_matrices_ ;
      Selection           selection_ ;
  } ; // selection


} } // namespace CORK::coefficient_matrices

#endif
