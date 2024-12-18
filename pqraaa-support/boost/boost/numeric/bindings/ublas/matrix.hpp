//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_UBLAS_MATRIX_HPP
#define BOOST_NUMERIC_BINDINGS_UBLAS_MATRIX_HPP

#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/end.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
#include <boost/numeric/bindings/detail/offset.hpp>
#include <boost/numeric/bindings/ublas/detail/convert_to.hpp>
#include <boost/numeric/bindings/ublas/storage.hpp>
#include <boost/numeric/bindings/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename F, typename A, typename Id, typename Enable >
struct adaptor< ::boost::numeric::ublas::matrix< T, F, A >, Id, Enable > {

    typedef typename copy_const< Id, T >::type value_type;
    typedef typename convert_to< tag::data_order, F >::type data_order;
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
        return id.size1();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.size2();
    }

    static value_type* begin_value( Id& id ) {
        return bindings::begin_value( id.data() );
    }

    static value_type* end_value( Id& id ) {
        return bindings::end_value( id.data() );
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return id.size2();
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return id.size1();
    }

};

template< typename T, std::size_t M, std::size_t N, typename F, typename Id, typename Enable >
struct adaptor< ::boost::numeric::ublas::bounded_matrix< T, M, N, F >, Id, Enable > {

    typedef typename copy_const< Id, T >::type value_type;
    typedef typename convert_to< tag::data_order, F >::type data_order;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::linear_array >,
        std::pair< tag::data_order, data_order >,
        std::pair< tag::stride_type<1>,
            typename if_row_major< data_order, mpl::int_<N>, tag::contiguous >::type >,
        std::pair< tag::stride_type<2>,
            typename if_row_major< data_order, tag::contiguous, mpl::int_<M> >::type >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size1();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.size2();
    }

    static value_type* begin_value( Id& id ) {
        return bindings::begin_value( id.data() );
    }

    static value_type* end_value( Id& id ) {
        return bindings::end_value( id.data() );
    }

};

template< typename T, std::size_t M, std::size_t N, typename Id, typename Enable >
struct adaptor< ::boost::numeric::ublas::c_matrix< T, M, N >, Id, Enable > {

    typedef typename copy_const< Id, T >::type value_type;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::size_type<2>, std::ptrdiff_t >,
        std::pair< tag::matrix_type, tag::general >,
        std::pair< tag::data_structure, tag::linear_array >,
        std::pair< tag::data_order, tag::row_major >,
        std::pair< tag::stride_type<1>, mpl::int_<N> >,
        std::pair< tag::stride_type<2>, tag::contiguous >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size1();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.size2();
    }

    static value_type* begin_value( Id& id ) {
        return id.data();
    }

    static value_type* end_value( Id& id ) {
        return id.data() + offset( id, id.size1(), id.size2() );
    }

};


} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
