//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_DETAIL_POD_HPP
#define BOOST_NUMERIC_BINDINGS_DETAIL_POD_HPP

#include <boost/numeric/bindings/is_numeric.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/copy_const.hpp>
#include <utility>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename Id >
struct adaptor< T, Id, typename std::enable_if< is_numeric<T>::value >::type > {

    typedef typename copy_const< Id, T >::type value_type;
/*    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::scalar >,
        mpl::pair< tag::size_type<1>, mpl::int_<1> >,
        mpl::pair< tag::data_structure, tag::linear_array >
    > property_map;
    */
    typedef std::tuple< std::pair< tag::value_type, value_type >
                      , std::pair< tag::entity, tag::scalar >
                      , std::pair< tag::size_type<1>, std::integral_constant<int,1> >
                      , std::pair< tag::data_structure, tag::linear_array >
                      > property_map ;

    static value_type* begin_value( Id& t ) {
        return &t;
    }

};

template< typename T, std::size_t N, typename Id >
struct adaptor< T[N], Id, typename std::enable_if< is_numeric<T>::value >::type > {

    typedef typename copy_const< Id, T >::type value_type;
    /*
    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::vector >,
        mpl::pair< tag::size_type<1>, mpl::int_<N> >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::stride_type<1>, tag::contiguous >
    > property_map;
    */
    typedef std::tuple< std::pair< tag::value_type, value_type >
                      , std::pair< tag::entity, tag::vector >
                      , std::pair< tag::size_type<1>, std::integral_constant<int,N> >
                      , std::pair< tag::data_structure, tag::linear_array >
                      , std::pair< tag::stride_type<1>, tag::contiguous >
                      > property_map ;

    static value_type* begin_value( Id& t ) {
        return &t[0];
    }

    static value_type* end_value( Id& t ) {
        return &t[N];
    }

};

template< typename T, std::size_t M, std::size_t N, typename Id >
struct adaptor< T[M][N], Id, typename std::enable_if< is_numeric<T>::value >::type > {

    typedef typename copy_const< Id, T >::type value_type;
    /*
    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::matrix >,
        mpl::pair< tag::size_type<1>, mpl::int_<M> >,
        mpl::pair< tag::size_type<2>, mpl::int_<N> >,
        mpl::pair< tag::matrix_type, tag::general >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::data_order, tag::row_major >,
        mpl::pair< tag::stride_type<1>, mpl::int_<N> >,
        mpl::pair< tag::stride_type<2>, tag::contiguous >
    > property_map;
    */
    typedef std::tuple< std::pair< tag::value_type, value_type >
                      , std::pair< tag::entity, tag::matrix >
                      , std::pair< tag::size_type<1>, std::integral_constant<int,M> >
                      , std::pair< tag::size_type<2>, std::integral_constant<int,N> >
                      , std::pair< tag::matrix_type, tag::general >
                      , std::pair< tag::data_structure, tag::linear_array >
                      , std::pair< tag::data_order, tag::row_major >
                      , std::pair< tag::stride_type<1>, std::integral_constant<int,N> >
                      , std::pair< tag::stride_type<2>, tag::contiguous >
                      > property_map ;

    static value_type* begin_value( Id& t ) {
        return &t[0][0];
    }

    static value_type* end_value( Id& t ) {
        return &t[M][N];
    }

};

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
