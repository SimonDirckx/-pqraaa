#ifndef glas2_mumps_sparse_hpp
#define glas2_mumps_sparse_hpp

#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/convert_to.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <type_traits>

namespace glas2 {
  namespace bindings_detail {

    template <typename E, typename EnableIf=void>
    struct inherits_from_contiguous_matrix
    : std::false_type
    {};

    template <typename E>
    struct inherits_from_contiguous_matrix< E, typename std::enable_if< is<DenseMatrix,E>::value >::type >
    : std::is_base_of< contiguous_matrix<typename E::value_type,typename E::size_type, typename E::orientation>, E >
    {};

  }
}

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template<>
struct convert_to< bindings::tag::data_order, glas2::row_major > {
  typedef bindings::tag::row_major type;
};

template<>
struct convert_to< bindings::tag::data_order, glas2::column_major > {
  typedef bindings::tag::column_major type;
};

template< typename E, typename Id >
struct adaptor< E, Id
              , typename std::enable_if< glas2::bindings_detail::inherits_from_contiguous_matrix<E>::value> ::type
              > {

    typedef typename copy_const< Id, typename E::value_type >::type               value_type;
    typedef typename convert_to< tag::data_order, typename E::orientation >::type data_order;
    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::matrix >,
        mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
        mpl::pair< tag::size_type<2>, std::ptrdiff_t >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::data_order, data_order >,
        mpl::pair< tag::stride_type<1>,
            typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
        mpl::pair< tag::stride_type<2>,
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


} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
