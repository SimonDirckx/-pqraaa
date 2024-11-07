//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_type_transformation_hpp
#define glas2_matrix_type_transformation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/type/transformation.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename M, typename Op>
  class transformation< M, Op
                        , typename std::enable_if< is<DenseMatrix,M>::value >::type
                        >
  {
    public:
      explicit transformation( M const& m )
      : matrix_( m )
      {}

    public:
      typedef typename M::size_type                         size_type ;
      typedef decltype( Op() ( typename M::value_type() ) ) value_type ;

      size_type num_rows() const {
        return matrix_.num_rows() ;
      }

      size_type num_columns() const {
        return matrix_.num_columns() ;
      }

      value_type operator() ( size_type i, size_type j ) const { return op_( matrix_(i,j) ) ; }
      value_type operator() ( size_type i, size_type j ) { return op_( matrix_(i,j) ) ; } // ?? REFERENCE

    public:
      template <typename I1, typename I2>
      typename std::enable_if< !(std::is_integral<I1>::value && std::is_integral<I2>::value)
                             , typename matrix_selection< transformation, I1, I2 >::result_type
                             >::type operator()( I1 const& s1, I2 const& s2 ) const {
        return matrix_selection< transformation, I1, I2 >::apply( *this, s1, s2 ) ;
      }
      
    public:
      transformation& operator=( transformation const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

      template <typename E>
      transformation& operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

    public:
      M matrix() const { return matrix_; }

    private:
      M  matrix_ ;
      Op op_ ;
  } ;

  template <typename M, typename Op, typename S1, typename S2>
  struct matrix_selection< transformation<M,Op>, S1, S2, typename std::enable_if< !(std::is_integral<S1>::value && std::is_integral<S2>::value) >::type > {
    typedef typename matrix_selection<M,S1,S2>::result_type temp_type ;
    typedef transformation< temp_type, Op >                 result_type ;

    static result_type apply( transformation<M,Op> m, S1 const& s1, S2 const& s2 ) {
      return result_type( m.matrix()( s1, s2 ) ) ;
    }
  } ;


  template <typename M, typename Op>
  struct glas_concept< transformation<M,Op>
                , typename std::enable_if< is<DenseMatrix,M>::value>::type
                >
  : DenseMatrix
  {} ;

} // namespace glas2

#endif
