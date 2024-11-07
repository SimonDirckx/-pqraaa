//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_view_sparse_vector_view_hpp
#define glas2_sparse_view_sparse_vector_view_hpp

#include <glas2/sparse/concept/sparse_vector.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename Data, typename Indices, int IndexBase>
  class sparse_vector_view {
    public:
      typedef typename Data::value_type    value_type ;
      typedef typename Indices::value_type index_type ;
      typedef typename Indices::size_type  size_type ;

      static index_type const index_base = IndexBase ;
      typedef Data     data_type ;
      typedef Indices  indices_type ;

    public:
      sparse_vector_view( data_type data, indices_type indices, index_type size )
      : data_( data )
      , indices_( indices )
      , size_( size )
      {}

      sparse_vector_view( sparse_vector_view const& that )
      : data_( that.data_ )
      , indices_( that.indices_ )
      , size_( that.size_ )
      {}

      // Vector
      size_type size() const { return size_ ; }
      size_type num_nz() const {
        assert( data_.size()==indices_.size() ) ;
        return data_.size() ;
      }

      sparse_vector_view& operator=( sparse_vector_view const& that ) {
        assert( that.size()==size() ) ;
        assert( that.num_nz()==num_nz() ) ;
        data_ = that.data_ ;
        indices_ = that.indices_ ;
        return *this ;
      }

      template <typename E>
      sparse_vector_view& operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

      data_type data() { return data_ ; }
      data_type data() const { return data_ ; }
      indices_type index() const { return indices_ ; }

    private:
      data_type    data_ ;
      indices_type indices_ ;
      index_type   size_ ;
  } ;


  template <typename Data, typename Indices, int IndexBase>
  struct glas_concept< sparse_vector_view< Data, Indices, IndexBase > >
  : SparseVector
  {};

} // namespace glas2


#endif
