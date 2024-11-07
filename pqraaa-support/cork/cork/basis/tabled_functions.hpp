//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_tabled_functions_hpp
#define cork_basis_tabled_functions_hpp

#include <cork/utility/value_type_for.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/basis/iterator.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace basis {

  template <typename X, typename Y>
  class tabled_functions
  {
    private:
      static_assert( glas2::is< glas2::DenseVector, X >::value ) ;
      static_assert( glas2::is< glas2::DenseMatrix, Y >::value ) ;

    public:
      typedef typename std::decay< X >::type x_type ;
      typedef typename std::decay< Y >::type y_type ;
      typedef typename x_type::size_type size_type ;

      template <typename T>
      using value_type_for = typename y_type::value_type ;

      template <typename T>
      using has_value_type_for = std::is_convertible<T, typename x_type::value_type> ;

    public:
      explicit tabled_functions( x_type const& x, y_type const& y )
      : x_( x )
      , y_( y )
      {
        assert( x_.size() == y_.num_rows() ) ;
      }

    public:
      size_type num_terms() const {
        return y_.num_columns() ;
      }

    public:
      template <typename ValueType, typename FunctionValues>
      void evaluate( ValueType const& arg, FunctionValues values ) const {
        auto i = glas2::max_ind( -abs_squared(x_-arg) ) ;
        values = y_( i, glas2::all() ) ;
      } // evaluate

    public:
      x_type const& x() const { return x_ ; }
      y_type const& y() const { return y_ ; }

    private:
      X x_ ;
      Y y_ ;
  } ; // class tabled_functions

  template <typename T>
  struct is_tabled_functions
  : std::false_type
  {} ;

  template <typename X, typename YTable>
  struct is_tabled_functions< tabled_functions<X, YTable> >
  : std::true_type
  {} ;

  template <typename X, typename YTable>
  decltype(auto) make_tabled_functions( X const& x, YTable const& y ) { return tabled_functions<X, YTable>( x, y ) ; }

} } // namespace CORK::basis


namespace CORK {

  template <typename X, typename Y>
  struct value_type< basis::tabled_functions<X,Y> >
  : value_type< typename std::decay< Y >::type >
  {} ;

  template <typename T, typename X, typename Y>
  struct value_type_for< T, basis::tabled_functions<X,Y> >
  : value_type< basis::tabled_functions<X,Y> >
  {} ;

} // namespace CORK

#endif
