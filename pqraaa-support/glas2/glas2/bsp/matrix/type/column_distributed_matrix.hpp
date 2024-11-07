//  (C) Copyright Karl Meerbergen & Albert-Jan Yzelman 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)
//  Modified by Zifan Liu 2015-04

#ifndef glas2_bsp_matrix_type_column_distributed_matrix_hpp
#define glas2_bsp_matrix_type_column_distributed_matrix_hpp

#include <glas2/concept/concept.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/matrix/view/matrix_selection.hpp>
#include <glas2/backend/current_backend.hpp>
#include <glas2/bsp/backend/default/matrix/assign.hpp>
#include <glas2/bsp/matrix/concept/bsp_column_matrix.hpp>
#include <glas2/bsp/vector/type/distributed_vector.hpp>
#include <type_traits>
#include <cassert>


namespace glas2 { namespace bsp {

  template <typename M, typename D>
  class column_distributed_matrix {
    public:
      typedef typename M::value_type value_type ;
      typedef typename M::size_type  size_type ;

    public:
      column_distributed_matrix( M const& m, D const& d )
      : m_( m )
      , d_( d )
      {}

    public:

      //Copy reference
      column_distributed_matrix( column_distributed_matrix const& that )
	      : m_(that.m_)
	      , d_( that.d_ )
	  {}

    public:
      typedef M local_type ;
      typedef D distribution_type ;
      local_type  local() const& { return m_ ; }
      distribution_type const& distribution() const { return d_ ; }

      size_type num_columns() const { return m_.num_columns() ; }
      size_type num_rows() const {return m_.num_rows();}
      /*
    public:
      column_distributed_matrix& operator=( column_distributed_matrix const& that ) {
        assert( d_==that.d_ ) ;
        m_ = that.m_ ;
        return *this ;
      }

      template <typename E>
      column_distributed_matrix& operator=( E const& that ) {
        assign( current_backend(), *this, that ) ;
        return *this ;
      }
*/
      //size_type num_rows() const { return m_.num_rows() ; }

    public: // We can select rows for this matrix
      template <typename I1>
      typename glas2::matrix_selection< bsp::column_distributed_matrix<M,D>, I1, glas2::all >::result_type operator()( I1 const& s1, glas2::all ) {
        return glas2::matrix_selection< bsp::column_distributed_matrix<M,D>, I1, glas2::all >::apply( *this, s1, glas2::all() ) ;
      }

    private:
      M        m_ ;
      D const& d_ ;
  } ;

} } // namespace glas2::bsp

namespace glas2 {

  template <typename M, typename D>
  struct concept< bsp::column_distributed_matrix<M,D> >
  : bsp::BSPColumnMatrix
  {} ;

  template <typename M, typename D, typename I1>
  struct matrix_selection< bsp::column_distributed_matrix<M,D>, I1, glas2::all
                         , typename std::enable_if< std::is_integral<I1>::value >::type > {
    typedef bsp::distributed_vector< typename matrix_selection< M, I1, glas2::all >::result_type, D > result_type ;

    static result_type apply( bsp::column_distributed_matrix<M,D> m, I1 i1, glas2::all ) {
      return result_type( m.local()(i1,glas2::all()), m.distribution() ) ;
    }
  } ;

  template <typename M, typename D, typename R1>
  struct matrix_selection< bsp::column_distributed_matrix<M,D>, R1, glas2::all
                         , typename std::enable_if< is<DenseVector,R1>::value >::type > {
    typedef bsp::column_distributed_matrix< typename matrix_selection< M, R1, glas2::all >::result_type, D > result_type ;

    static result_type apply( R1 const& r1, bsp::column_distributed_matrix<M,D> m, glas2::all ) {
      return result_type( m.local()(r1,glas2::all()), m.distribution() ) ;
    }
  } ;

} // namespace glas

#endif
