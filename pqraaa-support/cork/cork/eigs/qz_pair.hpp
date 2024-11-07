//  (C) Copyright Karl Meerbergen 2023.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_qz_pair_hpp
#define cork_eigs_qz_pair_hpp

#include <glas2/matrix.hpp>

namespace CORK { namespace eigs {

  template <typename T>
  class qz_pair {
    public:
      typedef T                                value_type ;
      typedef glas2::shared_matrix<value_type> matrix_type ;
      typedef typename matrix_type::size_type  size_type ;

    public:
      qz_pair( size_type n )
      : A_( n, n )
      , B_( n, n )
      {}

    public:
      auto size() const {
        return A_.num_rows() ;
      }

    public:
      matrix_type A_ ;
      matrix_type B_ ;
  } ;

} } // namespace CORK::eigs

#endif
