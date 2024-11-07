#ifndef glas2_bindings_matrix_adaptor_hpp
#define glas2_bindings_matrix_adaptor_hpp

#include <glas2/bindings/detail/data_order.hpp>
#include <glas2/matrix/concept/contiguous_dense_matrix.hpp>
#include <glas2/matrix/concept/strided_dense_matrix.hpp>
#include <glas2/matrix/concept/double_strided_dense_matrix.hpp>
#include <glas2/matrix/algorithm/upper.hpp>
#include <glas2/matrix/algorithm/lower.hpp>
#include <glas2/concept/is.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename E, typename Id >
struct adaptor< E, Id
              , typename std::enable_if< glas2::is<glas2::ContiguousDenseMatrix,E>::value> ::type
              > {

    typedef typename copy_const< Id, typename E::value_type >::type               value_type;
    typedef typename convert_to< tag::data_order, typename E::orientation >::type data_order;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::linear_array >,
        std::pair< tag::data_order, data_order >,
        std::pair< tag::stride_type<1>,
            typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
        std::pair< tag::stride_type<2>,
            typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.num_rows();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.num_columns();
    }

    static value_type* begin_value( Id& id ) {
        return id.ptr();
    }

    static value_type* end_value( Id& id ) {
        return id.ptr() + id.num_rows() * id.num_columns();
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return id.num_columns();
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return id.num_rows();
    }

};


namespace detail {
  template <typename Id>
  std::ptrdiff_t stride1( tag::row_major, const Id& id ) { return id.stride(); }

  template <typename Id>
  std::ptrdiff_t stride2( tag::row_major, const Id& id ) { return id.num_rows(); }

  template <typename Id>
  std::ptrdiff_t stride1( tag::column_major, const Id& id ) { return id.num_columns(); }

  template <typename Id>
  std::ptrdiff_t stride2( tag::column_major, const Id& id ) { return id.stride(); }

}

template< typename E, typename Id >
struct adaptor< E, Id
              , typename std::enable_if< glas2::is<glas2::StridedDenseMatrix,E>::value> ::type
              > {

    typedef typename copy_const< Id, typename E::value_type >::type               value_type;
    typedef typename convert_to< tag::data_order, typename E::orientation >::type data_order;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::linear_array >,
        std::pair< tag::data_order, data_order >,
        std::pair< tag::stride_type<1>,
            typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
        std::pair< tag::stride_type<2>,
            typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.num_rows();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.num_columns();
    }

    static value_type* begin_value( Id& id ) {
        return id.ptr();
    }

    static value_type* end_value( Id& id ) {
        return id.ptr() + stride1(id)*stride2(id) ;
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return detail::stride1( data_order(), id ) ;
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return detail::stride2( data_order(), id ) ;
    }

};



template< typename E, typename Id >
struct adaptor< E, Id
              , typename std::enable_if< glas2::is<glas2::DoubleStridedDenseMatrix,E>::value> ::type
              > {

    typedef typename copy_const< Id, typename E::value_type >::type               value_type;
    typedef typename convert_to< tag::data_order, typename E::orientation >::type data_order;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::linear_array >,
        std::pair< tag::data_order, data_order >,
        std::pair< tag::stride_type<1>,
            typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
        std::pair< tag::stride_type<2>,
            typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.num_rows();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.num_columns();
    }

    static value_type* begin_value( Id& id ) {
        return id.ptr();
    }

    static value_type* end_value( Id& id ) {
        return id.ptr() + std::max( id.stride_rows(), id.num_rows() ) * std::max( id.stride_columns(), id.num_columns() ) ;
    }

    static std::ptrdiff_t stride1( const Id& id ) {
      assert( id.stride_rows()==1 || id.stride_columns()==1 ) ;
      return id.stride_rows() ;
    }

    static std::ptrdiff_t stride2( const Id& id ) {
      assert( id.stride_rows()==1 || id.stride_columns()==1 ) ;
      return id.stride_columns() ;
    }

};


template< typename E, typename Id >
struct adaptor< glas2::lower_view<E>, Id
              , typename std::enable_if< glas2::is<glas2::DenseMatrix,E>::value> ::type
              > {

    typedef typename copy_const< Id, typename E::value_type >::type               value_type;
    typedef typename convert_to< tag::data_order, typename E::orientation >::type data_order;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::triangular_array >,
        std::pair< tag::data_side, tag::lower >,
        std::pair< tag::data_order, data_order >,
        std::pair< tag::stride_type<1>,
            typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
        std::pair< tag::stride_type<2>,
            typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.num_rows();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.num_columns();
    }

    static value_type* begin_value( Id& id ) {
        return id.matrix().ptr();
    }

    static value_type* end_value( Id& id ) {
        return id.matrix().ptr() + boost::numeric::bindings::stride1(id)*boost::numeric::bindings::stride2(id) ;
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return boost::numeric::bindings::stride1(id.matrix()) ;
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return boost::numeric::bindings::stride2(id.matrix()) ;
    }

};



template< typename E, typename Id >
struct adaptor< glas2::upper_view<E>, Id
              , typename std::enable_if< glas2::is<glas2::DenseMatrix,E>::value> ::type
              > {

    typedef typename copy_const< Id, typename E::value_type >::type               value_type;
    typedef typename convert_to< tag::data_order, typename E::orientation >::type data_order;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::triangular_array >,
        std::pair< tag::data_side, tag::upper >,
        std::pair< tag::data_order, data_order >,
        std::pair< tag::stride_type<1>,
            typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
        std::pair< tag::stride_type<2>,
            typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.num_rows();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.num_columns();
    }

    static value_type* begin_value( Id& id ) {
        return id.matrix().ptr();
    }

    static value_type* end_value( Id& id ) {
        return id.matrix().ptr() + boost::numeric::bindings::stride1(id)*boost::numeric::bindings::stride2(id) ;
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return boost::numeric::bindings::stride1(id.matrix()) ;
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return boost::numeric::bindings::stride2(id.matrix()) ;
    }

};


} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
