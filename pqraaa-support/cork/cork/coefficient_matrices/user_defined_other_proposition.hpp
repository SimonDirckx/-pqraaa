//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_user_defined_hpp
#define cork_coefficient_matrices_user_defined_hpp

#include <cork/coefficient_matrices/coefficient_matrices.hpp>
#include <cork/coefficient_matrices/user_defined_matrix.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <cassert>
#include <type_traits>
#include <glas2/vector.hpp>

namespace CORK { namespace coefficient_matrices {

  template <typename SequenceMatrices>
  class user_defined
  {
    public:
      typedef typename std::decay< SequenceMatrices>::type                sequence_type ;
      typedef typename sequence_type::value_type                          matrix_type ;
      typedef typename matrix_type::value_type                            value_type ;
      typedef typename matrix_type::size_type                             size_type ;
      typedef size_type                                                   grade_type ;

    public:
      user_defined( SequenceMatrices const& sequence )
      : sequence_( sequence )
      {
        assert( num_matrices()!=0 ) ;
      }

    public:
      size_type num_rows() const { return sequence_[0].num_rows() ; }
      size_type num_columns() const { return sequence_[0].num_columns() ; }

      grade_type num_matrices() const { return sequence_.size() ; }

    public:
      SequenceMatrices const& coefficient_matrices() const { return sequence_ ; }

      auto operator()(grade_type i) const { return sequence_[i] ; }

    public:
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        static_assert( glas2::is<glas2::ContiguousDenseVector,W>::value, "W should be a contiguous vector" ) ;
        assert( x.size()==num_rows() ) ;
        assert( w.size()==num_rows() ) ;
        assert( i>=0 && i<num_matrices() ) ;
        sequence_[i].multiply_add( x, w ) ;
      }

    private:
      // We do not provide accumulate

    private:
      SequenceMatrices sequence_ ;
  } ; // user_defined


  template <typename SequenceMatrices>
  user_defined< SequenceMatrices > make_user_defined( SequenceMatrices const& sequence ) { return user_defined<SequenceMatrices>( sequence ) ; }

  template <typename SequenceMatrices>
  user_defined< SequenceMatrices const& > make_user_defined_lvalue_reference( SequenceMatrices const& sequence ) { return user_defined<SequenceMatrices const&>( sequence ) ; }

  template <typename SequenceMatrices>
  struct coefficient_matrices_traits< SequenceMatrices, typename std::enable_if< is_user_defined_matrix<typename SequenceMatrices::value_type>::value >::type > {
    typedef user_defined< SequenceMatrices > type ;

    static type apply( SequenceMatrices const& sequence ) {
      return type( sequence ) ;
    }
  } ;

} } // namespace CORK::coefficient_matrices

#endif
