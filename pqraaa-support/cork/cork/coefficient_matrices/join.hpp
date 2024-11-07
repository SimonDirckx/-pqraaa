//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_join_hpp
#define cork_coefficient_matrices_join_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  //
  // List of matrices of different types
  // CoefficidentMatrices are:
  //   [Front, Tail]
  // where Front, Tail are CoefficientMatrices
  //
  template <typename Front, typename Tail>
  class join
  {
    public:
      typedef typename std::decay<Front>::type                              front_type ;
      typedef typename std::decay<Tail>::type                               tail_type ;
      typedef typename std::common_type< typename front_type::value_type
                                       , typename tail_type::value_type
                                       >::type                              value_type ;
      typedef typename std::common_type< typename front_type::size_type
                                       , typename tail_type::size_type
                                       >::type                              size_type ;
      typedef typename tail_type::size_type                                 grade_type ;

    public:
      join( Front front, Tail tail )
      : front_( front )
      , tail_( tail )
      {
        assert( front_.num_rows() == tail_.num_rows() ) ;
        assert( front_.num_columns() == tail_.num_columns() ) ;
      }

    public:
      size_type num_rows() const { return front_.num_rows() ; }
      size_type num_columns() const { return front_.num_columns() ; }

      grade_type num_matrices() const { return front_.num_matrices() + tail_.num_matrices() ; }

    public:
      front_type const& front() const { return front_ ; }
      tail_type const& tail() const { return tail_ ; }

    public:
      template <typename M>
      void fill( int i, M m ) const {
        if (i<front_.num_matrices())
          front_.fill( i, m ) ;
        else
          tail_.fill( i-front_.num_matrices(), m ) ;
      } // fill()

    public:
      template <typename X, typename W>
      void multiply_add( size_type i, X const& x, W& w ) const {
        assert( i>=0 && i< num_matrices() ) ;
        if (i<front_.num_matrices()) 
          front_.multiply_add( i, x, w ) ;
        else
          tail_.multiply_add( i-front_.num_matrices(), x, w ) ;
      } // apply_scheduled()

      template <typename Coefs, typename Matrix>
      void accumulate( Coefs const& coefs, Matrix A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
        assert( coefs.size()==num_matrices() ) ;

        front_.accumulate( coefs(glas2::range(0,front_.num_matrices())), A ) ;
        tail_.accumulate( coefs(glas2::range_from_end(front_.num_matrices(),0)), A ) ;
      } // accumulate()

    private:
      Front  front_ ;
      Tail   tail_ ;
  } ; // join

  template <typename Front, typename Tail>
  join<Front const&, Tail const&> make_join_lvalue_reference( Front const& d, Tail const& t ) { return join<Front const&, Tail const&>(d,t); }

  template <typename Front, typename Tail>
  join<Front, Tail> make_join( Front const& d, Tail const& t ) { return join<Front, Tail>(d,t); }

} } // namespace CORK::coefficient_matrices

#endif
