//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_matrix_assign_hpp
#define glas2_backend_default_backend_matrix_assign_hpp

#include <glas2/backend/default_backend/default_backend.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>
#include <iostream>

namespace glas2 {

/*  namespace detail_assign {
    template <typename To, typename From>
    To assign_dense_matrix( std::true_type, To to, From const& from ) {
      assert( from.num_rows() == to.num_rows() ) ;
      assert( from.num_columns() == to.num_columns() ) ;

      auto it_to = to.iterate() ;
      auto it_from = from.iterate() ;
      for (typename To::size_type i=0; i<it_to.size_1(); ++i) {
        for (typename To::size_type j=0; j<it_to.size_2(); ++j) {
          to[*it_to] = from[*it_from] ;
          it_from.index_2_pp() ;
          it_to.index_2_pp() ;
        }
        it_from.index_1_pp() ;
        it_to.index_1_pp() ;
      }
      return to ;
    }

    template <typename To, typename From>
    To assign_dense_matrix( std::false_type, To to, From const& from ) {
      std::cerr << "Assignment of matrices with different orientations is not implemented" << std::endl ;
      assert( false ) ;
    }
  } // namespace detail_assign

  template <typename To, typename From>
  typename std::enable_if< boost::mpl::and_< is<DenseMatrix,To>, is<DenseMatrix,From> >::value, To >::type assign( To to, From const& from ) {
    typedef std::is_same< typename From::orientation, typename To::orientation > same_orientation ;
    detail_assign::assign_dense_matrix( same_orientation(), to, from ) ;
    return to ;
  }
  */

  template <typename To, typename From>
  typename std::enable_if< is<DenseMatrix,To>::value && is<DenseMatrix,From>::value, To& >::type assign( default_backend, To& to, From const& from ) {
    assert( from.num_rows()==to.num_rows() );
    assert( from.num_columns()==to.num_columns() );

    for (typename To::size_type i=0; i<to.num_rows(); ++i) {
      for (typename To::size_type j=0; j<to.num_columns(); ++j) {
        to(i,j) = from(i,j) ;
      }
    }
    return to ;
  }

} // namespace glas2

#endif
