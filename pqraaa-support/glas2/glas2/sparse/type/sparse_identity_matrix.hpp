//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_sparse_type_sparse_identity_matrix_hpp
#define glas2_sparse_type_sparse_identity_matrix_hpp

#include <glas2/vector/type/range.hpp>
#include <glas2/vector/type/slice.hpp>
#include <glas2/vector/type/all.hpp>
#include <glas2/vector/type/constant_vector.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 {

  template <typename T=double, int IndexBase=0>
  class sparse_identity_matrix {
    public:
      typedef T         value_type ;
      typedef int       size_type ;
      typedef size_type index_type ;
      static const int index_base = IndexBase ;

    public:
      sparse_identity_matrix()
      : num_rows_(0)
      , num_columns_(0)
      {}

      sparse_identity_matrix( size_type m, size_type n )
      : num_rows_( m )
      , num_columns_( n )
      {}

      // Copy reference !!
      sparse_identity_matrix( sparse_identity_matrix const& that )
      : num_rows_( that.num_rows_ )
      , num_columns_( that.num_columns_ )
      {}

    public:
      void reset( size_type m, size_type n ) {
        num_rows_ = m ;
        num_columns_ = n ;
      }

    public:
      // Matrix
      size_type num_rows() const { return num_rows_ ; }
      size_type num_columns() const { return num_columns_ ; }
      size_type num_nz() const { return std::min(num_rows_, num_columns_) ; }

      typedef glas2::constant_vector<size_type,value_type> data_type ;
      data_type data() const{ return data_type(num_nz(), 1) ; }

      typedef glas2::range rows_type ;
      glas2::range rows() const{ return glas2::range(IndexBase,num_nz()+IndexBase);}

      typedef glas2::range columns_type ;
      glas2::range columns() const{ return glas2::range(IndexBase,num_nz()+IndexBase);}

      size_type row(size_type i) const {
        assert( i>=0 && i<num_nz() ) ;
        return i ; 
      }

      size_type column(size_type i) const {
        assert( i>=0 && i<num_nz() ) ;
        return i ; 
      }

    private:
      size_type num_rows_ ;
      size_type num_columns_ ;
  } ;


  template <typename T, int IndexBase>
  struct glas_concept< sparse_identity_matrix<T,IndexBase> >
  : CoordinateSparseMatrix
  {};

} // namespace glas2


#endif
