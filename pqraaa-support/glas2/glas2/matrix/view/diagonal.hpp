#ifndef glas2_matrix_view_diagonal_view_hpp
#define glas2_matrix_view_diagonal_view_hpp

#include <glas2/concept/concept.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <cassert>

namespace glas2 {

  template <typename O>
  struct diagonal_view_orientation
  {} ;

  template <>
  struct diagonal_view_orientation< row_major >
  {
    typedef column_major type ;
  } ;

  template <>
  struct diagonal_view_orientation< column_major >
  {
    typedef row_major type ;
  } ;

  template <typename M>
  class diagonal_view {
    public:
      typedef typename M::value_type value_type ;
      typedef typename M::size_type  size_type ;
      typedef typename diagonal_view_orientation< typename M::orientation >::type orientation ;

    public:
      diagonal_view( M m )
      : m_(m)
      {}

    public:
      size_type num_rows() const { return m_.num_columns() ; }
      size_type num_columns() const { return m_.num_rows() ; }

      value_type operator() ( size_type i, size_type j ) const {
        assert( i>=0 && i<num_rows() ) ;
        assert( j>=0 && j<num_columns() ) ;
        return m_(j,i) ;
      }
      value_type& operator() ( size_type i, size_type j ) {
        return m_(j,i) ;
      }

    /*public:
      typedef typename M::iterator iterator ;
      iterator iterate() const { return m_.iterate() ; }

      value_type operator[] ( size_type i ) const{ return m_[i] ; }
      value_type& operator[] ( size_type i ) { return m_[i] ; }
     */

    private:
      M m_ ;
  } ;

  template <typename M>
  struct glas_concept< diagonal_view<M> >
  : glas_concept<M>
  {} ;

} // namespace glas2

#endif
