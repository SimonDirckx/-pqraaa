//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_list_hpp
#define cork_coefficient_matrices_list_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  //
  // List of matrices of different types
  // CoefficidentMatrices are:
  //   [Matrix, Tail]
  // where Tail is CoefficientMatrices
  //
  template <typename Matrix, typename Tail>
  class list
  {
    public:
      typedef typename std::decay<Matrix>::type                             first_type ;
      typedef typename std::decay<Tail>::type                               tail_type ;
      typedef typename std::common_type< typename first_type::value_type
                                       , typename tail_type::value_type
                                       >::type                              value_type ;
      typedef typename std::common_type< typename first_type::size_type
                                       , typename tail_type::size_type
                                       >::type                              size_type ;
      typedef typename tail_type::size_type                                 grade_type ;

    public:
      list( Matrix first, Tail tail )
      : first_( first )
      , tail_( tail )
      {
        assert( first_.num_rows() == tail_.num_rows() ) ;
        assert( first_.num_ciolumns() == tail_.num_ciolumns() ) ;
      }

    public:
      size_type num_rows() const { return first_.num_rows() ; }
      size_type num_columns() const { return first_.num_columns() ; }

      grade_type num_matrices() const { return tail_.num_matrices()+1 ; }

    public:
      template <typename M>
      void fill( int i, M m ) const {
        if (i==0) {
          m = first_ ;
        } else {
          tail_.fill( i-1, m ) ;
        }
      }

    public:
      template <typename X, typename W>
      void multiply_add( size_type i, X const& x, W& w ) const {
        assert( i>=0 && i< num_matrices() ) ;
        if (i==0) 
          w += glas2::multiply( first_, x ) ;
        else
          tail_.multiply_add( i-1, x, w ) ;
      } // apply_scheduled()

      template <typename Coefs, typename Matrix>
      void accumulate( Coefs const& coefs, Matrix A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
        assert( coefs.size()==num_matrices() ) ;

        A += coefs(0) * first_ ;
        tail_.accumulate( coefs(glas2::range_from_end(1,0)), A ) ;
      } // accumulate()

    private:
      Matrix first_ ;
      Tail   tail_ ;
  } ; // list

  template <typename Matrix, typename Tail>
  list<Matrix const&, Tail const&> const& make_list_lvalue_reference( Matrix const& d, Tail const& t ) { return list<Matrices const&, Tail const&>(d,t); }

  template <typename Matrix, typename Tail>
  list<Matrix, Tail> const& make_list( Matrix const& d, Tail const& t ) { return list<Matrices, Tail>(d,t); }

} } // namespace CORK::coefficient_matrices

#endif
