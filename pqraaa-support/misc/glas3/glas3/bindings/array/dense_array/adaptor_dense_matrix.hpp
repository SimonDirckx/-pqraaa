//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_bindings_array_dense_array_adaptor_dense_matrix_hpp
#define glas3_bindings_array_dense_array_adaptor_dense_matrix_hpp

#include <glas3/array/dense_array/concept/contiguous_dense_array.hpp>
//#include <glas2/matrix/algorithm/upper.hpp>
#include <glas3/concept/is.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/convert_to.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

//template<>
//struct convert_to< bindings::tag::data_order, glas2::row_major > {
//  typedef bindings::tag::row_major type;
//};
//
//template<>
//struct convert_to< bindings::tag::data_order, glas2::column_major > {
//  typedef bindings::tag::column_major type;
//};

template< typename E, typename Id >
struct adaptor< E, Id, typename std::enable_if< glas3::is<glas3::ContiguousDenseMatrix,E>::value> ::type > {

    typedef typename copy_const< Id, typename E::value_type >::type     value_type;
    typedef typename bindings::tag::column_major                        data_order;
    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::matrix >,
        mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
        mpl::pair< tag::size_type<2>, std::ptrdiff_t >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::data_order, data_order >,
        mpl::pair< tag::stride_type<1>, tag::contiguous >,
        mpl::pair< tag::stride_type<2>, std::ptrdiff_t >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.shape()[0] ;
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.shape()[1] ;
    }

    static value_type* begin_value( Id& id ) {
        return id.data_ptr().get() ;
    }

    static value_type* end_value( Id& id ) {
        return id.data_ptr().get() + id.size() ;
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return id.shape()[1] ;
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return id.shape()[0] ;
    }

};

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
