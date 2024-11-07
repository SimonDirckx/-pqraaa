#ifndef glas2_bindings_matrix_expression_hpp
#define glas2_bindings_matrix_expression_hpp

#include <glas2/bindings/detail/data_order.hpp>
#include <glas2/matrix/expression/unary_operation.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/concept/conjugate.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/conj.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename E, typename Id >
struct adaptor< glas2::unary_operation< E, glas2::conjugate >, Id > {

    typedef typename copy_const< Id, typename glas2::unary_operation< E, glas2::conjugate >::value_type >::type               value_type;
    typedef typename convert_to< tag::data_order, typename E::orientation >::type data_order;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::linear_array >,
        std::pair< tag::data_order, data_order >,
        typename std::conditional<
            is_same_at< E, tag::value_transform, tag::conjugate >::value,
            std::pair< tag::value_transform, std::false_type >,
            std::pair< tag::value_transform, tag::conjugate >
        >::type,
        std::pair< tag::stride_type<1>,
            typename if_row_major< data_order, std::ptrdiff_t, tag::contiguous >::type >,
        std::pair< tag::stride_type<2>,
            typename if_row_major< data_order, tag::contiguous, std::ptrdiff_t >::type >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return ::boost::numeric::bindings::size1(id.matrix());
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return ::boost::numeric::bindings::size2(id.matrix());
    }

    static value_type* begin_value( Id& id ) {
      return ::boost::numeric::bindings::begin_value(id.matrix());
    }

    static value_type* end_value( Id& id ) {
      return ::boost::numeric::bindings::end_value(id.matrix()) ;
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return ::boost::numeric::bindings::stride1(id.matrix());
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return ::boost::numeric::bindings::stride2(id.matrix());
    }

};




} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
