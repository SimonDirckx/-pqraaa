//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_algorithm_iostream_hpp
#define glas2_matrix_algorithm_iostream_hpp

#include <glas2/matrix/concept/matrix.hpp>
#include <glas2/type/nice_print.hpp>
#include <glas2/concept/is.hpp>
#include <iostream>
#include <type_traits>

namespace glas2 {

  template <typename V>
  typename std::enable_if< is<Matrix,V>::value, std::ostream&>::type operator<< ( std::ostream& s, V const& v ) {
    s << "(" << v.num_rows() << "," << v.num_columns() << ")[" ;
    for ( typename V::size_type i=0; i<v.num_rows(); ++i) {
      if (v.num_columns()>0) s << "[" << v(i,0) ;
      for (typename V::size_type j=1; j<v.num_columns(); ++j) { s << "," << v(i,j) ; }
      s << "];" ;
    }
    s << "]" ;
    return s ;
  }

  template <typename M>
  typename std::enable_if< is<Matrix,M>::value, std::ostream&>::type operator<< ( std::ostream& s, nice_print_type<M> const& m ) {
    s << "\n" ;
    for ( typename M::size_type i=0; i<m.e_.num_rows(); ++i) {
      for (typename M::size_type j=0; j<m.e_.num_columns(); ++j) { s << " " << m.e_(i,j) ; }
      s << "\n" ;
    }
    return s ;
  }

} // namespace glas2

#endif
