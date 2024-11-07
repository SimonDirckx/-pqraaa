//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_view_sparse_vector_view_void_hpp
#define glas2_sparse_view_sparse_vector_view_void_hpp

#include <glas2/sparse/view/sparse_vector_view.hpp>
#include <glas2/sparse/type/no_value_type_array.hpp>

namespace glas2 {

  template <typename Indices, int IndexBase>
  class sparse_vector_view< no_value_type_array, Indices, IndexBase > {
    public:
      typedef no_value_type                value_type ;
      typedef typename Indices::value_type index_type ;
      typedef typename Indices::size_type  size_type ;

      static index_type const index_base = IndexBase ;
      typedef Indices             indices_type ;
      typedef no_value_type_array data_type ;

    public:
      sparse_vector_view( indices_type indices, index_type size )
      : indices_( indices )
      , size_( size )
      {}

      sparse_vector_view( sparse_vector_view const& that )
      : indices_( that.indices_ )
      , size_( that.size_ )
      {}

      // Vector
      size_type size() const { return size_ ; }
      size_type num_nz() const {
        return indices_.size() ;
      }

      sparse_vector_view& operator=( sparse_vector_view const& that ) {
        assert( that.size()==size() ) ;
        assert( that.num_nz()==num_nz() ) ;
        indices_ = that.indices_ ;
        return *this ;
      }

      template <typename E>
      sparse_vector_view& operator=( E const& that ) {
        assign( *this, that ) ;
        return *this ;
      }

      data_type data() const { return no_value_type_array() ; }
      indices_type index() const { return indices_ ; }

    private:
      indices_type indices_ ;
      index_type   size_ ;
  } ;


  template <typename Indices, int IndexBase>
  class glas_concept< sparse_vector_view< no_value_type, Indices, IndexBase > >
  : SparseVector
  {};

} // namespace glas2


#endif
