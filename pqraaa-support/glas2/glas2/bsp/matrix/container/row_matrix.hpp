//  (C) Copyright Karl Meerbergen & Albert-Jan Yzelman 2014
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)
//  Modified by Zifan Liu

#ifndef glas2_bsp_matrix_container_row_matrix_hpp
#define glas2_bsp_matrix_container_row_matrix_hpp

#include <glas2/bsp/matrix/type/row_distributed_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_row_matrix.hpp>
#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 { namespace bsp {

  template <typename Distribution, typename T, typename O=column_major, typename S=std::ptrdiff_t>
  class row_matrix
  : public row_distributed_matrix< contiguous_matrix<T,S,O>, Distribution >
  {
    public:
      typedef row_distributed_matrix< contiguous_matrix<T,S,O>, Distribution > base_type ;
      typedef typename base_type::local_type                                   local_type ;
      typedef typename local_type::size_type                                   size_type ;
    public:
      row_matrix( size_type m, size_type n, Distribution const& distribution=Distribution() )
      : base_type( contiguous_matrix<T,S,O>( new T[m*distribution.getSize()], m, distribution.getSize() ), distribution )
      , owner_( true ) 
      {}

      // Copy reference !!
      row_matrix( row_matrix const& that )
      : base_type( that )
      , owner_( false )
      {}

      ~row_matrix() {
#ifndef NDEBUG
        if (owner_) std::cout << "Calling destructor BSP::ROW_MATRIX" << std::endl ;
#endif
        if (owner_) delete [] this->local().ptr() ;
      }
/*
    public:
      row_matrix operator=( row_matrix const& that ) {
        base_type(*this) = that ;
        return *this ;
      }

      template <typename E>
      row_matrix operator=( E const& that ) {
        base_type(*this) = that ;
        return *this ;
      }
      */
    private:
      bool owner_ ;
  } ;

} } // namespace glas2::bsp


namespace glas2 {

  template <typename D, typename T, typename O>
  struct concept< bsp::row_matrix<D,T,O> >
  : bsp::BSPRowMatrix
  {};

/*  template <typename Tag, typename T, typename O>
  struct matrix_transform< Tag, shared_matrix<T,O> >
  : matrix_transform< Tag, typename shared_matrix<T,O>::base_type >
  {} ;*/

} // namespace glas2

#endif
