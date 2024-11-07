#ifndef glas2_bindings_matrix_bind_matrix_hpp
#define glas2_bindings_matrix_bind_matrix_hpp

#include <glas2/matrix/concept/contiguous_dense_matrix.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/matrix/type/strided_matrix.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/num_rows.hpp>
#include <boost/numeric/bindings/num_columns.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_row_major.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename EnableIf=void>
  struct bind_matrix_selector {
    typedef typename boost::numeric::bindings::value_type<M>::type                           value_type ;
    typedef typename boost::numeric::bindings::result_of::size< M >::type                    size_type ;
    typedef typename std::conditional< boost::numeric::bindings::is_column_major<M>::value
                                     , column_major, row_major
                                     >::type                                                 orientation_type ;
    typedef strided_matrix< value_type, size_type, orientation_type >                        type ;

    static type apply( M& m ) {
      return type( boost::numeric::bindings::begin_value(m), boost::numeric::bindings::stride_major(m)
                 , boost::numeric::bindings::num_rows(m), boost::numeric::bindings::num_columns(m)
                 ) ;
    }
  } ;

  template <typename M>
  typename bind_matrix_selector< M >::type bind_matrix( M& m ) {
    return bind_matrix_selector< M >::apply( m ) ;
  }

} // namespace glas2

#endif
