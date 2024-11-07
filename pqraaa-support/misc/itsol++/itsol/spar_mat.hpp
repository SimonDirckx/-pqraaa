#ifndef ITSOL_SPAR_MAT_HPP
#define ITSOL_SPAR_MAT_HPP

#include <boost/numeric/bindings/index_base.hpp>
#include <boost/numeric/bindings/num_rows.hpp>
#include <boost/numeric/bindings/num_columns.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/end.hpp>
#include <boost/numeric/bindings/data_order.hpp>
#include <type_traits>
#include <vector>
#include <cassert>

namespace itsol {

  template <typename T>
  struct spar_mat_type {
    int   n_ ;
    std::vector<int>  nzcount_ ;
    std::vector<int*> ja_ ;
    std::vector<T*>   ma_ ;
  } ;

  template <typename M, typename T>
  typename std::enable_if< std::is_same< boost::numeric::bindings::tag::compressed_sparse
                                       , typename boost::numeric::bindings::detail::property_at< M, boost::numeric::bindings::tag::data_structure >::type
                                       >::value
                         >::type spar_mat( M& matrix, spar_mat_type<T>& sp_mat ) {
     typedef typename ::boost::numeric::bindings::result_of::index_base<M>::type index_base ;
     static_assert( index_base::value==0, "index_base should be 0" ) ;
     static_assert( boost::is_same<typename boost::numeric::bindings::value_type<M>::type,T>::value, "value_type of matrices should be the same" ) ;
     static_assert( boost::is_same<typename boost::numeric::bindings::result_of::data_order<M>::type,boost::numeric::bindings::tag::row_major>::value, "matrix should be row major oriented" ) ;

     //int n = boost::numeric::bindings::num_rows( matrix ) ;
     sp_mat.n_ = boost::numeric::bindings::detail::adaptor<M,M>::size1( matrix ) ;
     assert( (sp_mat.n_==::boost::numeric::bindings::detail::adaptor<M,M>::size2( matrix )) ) ;
     //assert( sp_mat.n_==boost::numeric::bindings::num_columns( matrix ) ) ;

     auto rows = boost::numeric::bindings::begin_compressed_index_major(matrix) ;
     sp_mat.nzcount_.resize( sp_mat.n_ ) ;
     sp_mat.ja_.resize( sp_mat.n_ ) ;
     sp_mat.ma_.resize( sp_mat.n_ ) ;
     for (int i=0; i<sp_mat.n_; ++i) {
       sp_mat.nzcount_[i] = rows[i+1] - rows[i] ;
       sp_mat.ja_[i] = boost::numeric::bindings::begin_index_minor(matrix) + rows[i] ;
       sp_mat.ma_[i] = boost::numeric::bindings::begin_value(matrix) + rows[i] ;
     }
  } // spar_mat()

  //
  // Coordinate matrix M should be sorted row wise!
  // This is not tested.
  //
  template <typename M, typename T>
  typename std::enable_if< std::is_same< boost::numeric::bindings::tag::coordinate_sparse
                                       , typename boost::numeric::bindings::detail::property_at< M, boost::numeric::bindings::tag::data_structure >::type
                                       >::value
                         >::type spar_mat( M& matrix, spar_mat_type<T>& sp_mat ) {
     typedef typename ::boost::numeric::bindings::result_of::index_base<M>::type index_base ;
     static_assert( index_base::value==0, "index_base should be 0" ) ;
     static_assert( boost::is_same<typename boost::numeric::bindings::value_type<M>::type,T>::value, "value_type of matrices should be the same" ) ;

     //int n = boost::numeric::bindings::num_rows( matrix ) ;
     sp_mat.n_ = boost::numeric::bindings::detail::adaptor<M,M>::size1( matrix ) ;
     assert( (sp_mat.n_==::boost::numeric::bindings::detail::adaptor<M,M>::size2( matrix )) ) ;
     //assert( sp_mat.n_==boost::numeric::bindings::num_columns( matrix ) ) ;

     int nz = boost::numeric::bindings::end_value(matrix) - boost::numeric::bindings::begin_value(matrix) ;
     auto index_major = boost::numeric::bindings::begin_index_major(matrix) ;
     auto index_minor = boost::numeric::bindings::begin_index_minor(matrix) ;
     auto value = boost::numeric::bindings::begin_value(matrix) ;

     sp_mat.nzcount_.resize( sp_mat.n_ ) ;
     sp_mat.ja_.resize( sp_mat.n_ ) ;
     sp_mat.ma_.resize( sp_mat.n_ ) ;
     int row_i = 0 ; int pos_i = 0 ;
     int i ;
     for (i=0; i<nz; ++i) {
       int row = index_major[i];
       while (row!=row_i) {
         sp_mat.nzcount_[row_i] = i - pos_i ;
         sp_mat.ja_[row_i] = index_minor + pos_i ;
         sp_mat.ma_[row_i] = value + pos_i ;
         ++row_i ;
         pos_i = i ;
       }
     }
     for ( ; row_i<sp_mat.n_; ++row_i) {
       sp_mat.nzcount_[row_i] = nz - pos_i ;
       sp_mat.ja_[row_i] = index_minor + pos_i ;
       sp_mat.ma_[row_i] = value + pos_i ;
     }
  } // spar_mat()

} // namespace itsol

#endif
