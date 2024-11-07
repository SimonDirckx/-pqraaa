//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bsp_matrix_expression_binary_operation_hpp
#define glas2_bsp_matrix_expression_binary_operation_hpp

#include <glas2/concept/is.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/expression/binary_operation.hpp>
#include <glas2/scalar/concept/scalar.hpp>
#include <glas2/bsp/matrix/concept/bsp_row_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_column_matrix.hpp>
#include <type_traits>

namespace glas2 {

  template <typename S, typename M, typename Op>
  class binary_operation< S, M, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPMatrix,M>::value>::type
                        >
  {
    public:
      typedef binary_operation< S, typename M::local_type, Op > local_type ;
      typedef typename M::distribution_type                     distribution_type ;

    public:
      binary_operation( S const& s, M const& v )
      : local_( s, v.local() )
      , distribution_( v.distribution() )
      {}

    public:
      distribution_type const& distribution() const { return distribution_ ; }

      typedef typename M::size_type size_type ;
      typename std::enable_if< is<bsp::BSPRowMatrix,M>::value, size_type > num_columns() const {return local_.num_columns() ; }
      typename std::enable_if< is<bsp::BSPColumnMatrix,M>::value, size_type > num_rows() const {return local_.num_rows() ; }

      local_type const& local() const { return local_ ; }

    private:
      local_type        local_ ;
      distribution_type distribution_ ;
  } ;

  template <typename S, typename M, typename Op>
  struct concept< binary_operation<S,M,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPRowMatrix,M>::value>::type
                >
  : bsp::BSPRowMatrix
  {} ;

  template <typename S, typename M, typename Op>
  struct concept< binary_operation<S,M,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPColumnMatrix,M>::value>::type
                >
  : bsp::BSPColumnMatrix
  {} ;



  template <typename M, typename S, typename Op>
  class binary_operation< M, S, Op
                        , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPMatrix,M>::value>::type
                        >
  {
    public:
      typedef binary_operation< typename M::local_type, S, Op > local_type ;
      typedef typename M::distribution_type                     distribution_type ;

    public:
      binary_operation( M const& v, S const& s )
      : local_( v.local(), s )
      , distribution_( v.distribution() )
      {}

    public:
      distribution_type const& distribution() const { return distribution_ ; }

      local_type const& local() const { return local_ ; }

      typedef typename M::size_type size_type ;
      typename std::enable_if< is<bsp::BSPRowMatrix,M>::value, size_type > num_columns() const {return local_.num_columns() ; }
      typename std::enable_if< is<bsp::BSPColumnMatrix,M>::value, size_type > num_rows() const {return local_.num_rows() ; }

    private:
      local_type        local_ ;
      distribution_type distribution_ ;
  } ;

  template <typename M, typename S, typename Op>
  struct concept< binary_operation<M,S,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPColumnMatrix,M>::value>::type
                >
  : bsp::BSPColumnMatrix
  {} ;

  template <typename M, typename S, typename Op>
  struct concept< binary_operation<M,S,Op>
                , typename std::enable_if< is<Scalar,S>::value && is<bsp::BSPRowMatrix,M>::value>::type
                >
  : bsp::BSPRowMatrix
  {} ;


  template <typename M1, typename M2, typename Op>
  class binary_operation< M1, M2, Op
                        , typename std::enable_if< is<bsp::BSPMatrix,M1>::value && is<bsp::BSPMatrix,M2>::value>::type
                        >
  {
    public:
      typedef binary_operation< typename M1::local_type, typename M2::local_type, Op > local_type ;
      typedef typename M1::distribution_type                                           distribution_type ;

    public:
      binary_operation( M1 const& v1, M2 const& v2 )
      : local_( v1.local(), v2.local() )
      , distribution_( v1.distribution() )
      {
        assert( distribution_==v2.distribution() ) ;
      }

    public:
      distribution_type const& distribution() const { return distribution_ ; }

      local_type const& local() const { return local_ ; }

      typedef typename local_type::size_type size_type ;
      typename std::enable_if< is<bsp::BSPRowMatrix,M1>::value, size_type > num_columns() const {return local_.num_columns() ; }
      typename std::enable_if< is<bsp::BSPColumnMatrix,M1>::value, size_type > num_rows() const {return local_.num_rows() ; }

    private:
      local_type        local_ ;
      distribution_type distribution_ ;
  } ;

  template <typename M1, typename M2, typename Op>
  struct concept< binary_operation<M1,M2,Op>
                , typename std::enable_if< is<bsp::BSPRowMatrix,M1>::value && is<bsp::BSPRowMatrix,M2>::value >::type
                >
  : bsp::BSPRowMatrix
  {} ;

  template <typename M1, typename M2, typename Op>
  struct concept< binary_operation<M1,M2,Op>
                , typename std::enable_if< is<bsp::BSPColumnMatrix,M1>::value && is<bsp::BSPColumnMatrix,M2>::value >::type
                >
  : bsp::BSPColumnMatrix
  {} ;


} // namespace glas2

#endif
