//Zifan Liu, zifan.liu@cs.kuleuven.be, 2015-05-15

#ifndef glas2_bsp_matrix_container_col_matrix_hpp
#define glas2_bsp_matrix_container_col_matrix_hpp

#include <glas2/bsp/matrix/type/column_distributed_matrix.hpp>
#include <glas2/bsp/matrix/concept/bsp_column_matrix.hpp>
#include <glas2/matrix/type/contiguous_matrix.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 { namespace bsp {

  template <typename Distribution, typename T, typename O=row_major, typename S=std::ptrdiff_t>
  class col_matrix
  : public column_distributed_matrix< contiguous_matrix<T,S,O>, Distribution >
  {
    public:
      typedef column_distributed_matrix< contiguous_matrix<T,S,O>, Distribution > base_type ;
      typedef typename base_type::local_type                                   local_type ;
      typedef typename local_type::size_type                                   size_type ;
    public:
      col_matrix( size_type m, size_type n, Distribution const& distribution=Distribution() )
      : base_type( contiguous_matrix<T,S,O>( new T[distribution.getSize()*n], distribution.getSize(), n ), distribution )
      , owner_( true ) 
      {}

      // Copy reference !!
      col_matrix( col_matrix const& that )
      : base_type( that )
      , owner_( false )
      {}

      ~col_matrix() {
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
  struct concept< bsp::col_matrix<D,T,O> >
  : bsp::BSPColumnMatrix
  {};

/*  template <typename Tag, typename T, typename O>
  struct matrix_transform< Tag, shared_matrix<T,O> >
  : matrix_transform< Tag, typename shared_matrix<T,O>::base_type >
  {} ;*/

} // namespace glas2

#endif
