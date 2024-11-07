//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_TRANS_HPP
#define BOOST_NUMERIC_BINDINGS_TRANS_HPP

#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/end.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/rank.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/numeric/bindings/has_linear_array.hpp>
#include <type_traits>
#include <functional>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename Conj >
struct trans_wrapper: std::reference_wrapper<T> {
    trans_wrapper( T& t ): std::reference_wrapper<T>( t ) {}
};

//
// In case of linear storage
//
template< typename T, typename Conj, typename Id, typename Enable >
struct adaptor< trans_wrapper<T, Conj>, Id, Enable > {

    typedef typename property_map_of< T >::type prop_of_T;
    typedef typename property_insert< T,

        // upgrade to at least a matrix
        std::pair<
            tag::entity, 
            tag::tensor< std::max( tag::matrix::value, rank< T >::value ) >
        >,

        // size1 <-> size2
        std::pair< tag::size_type<1>, typename result_of::size2< T >::type >,
        std::pair< tag::size_type<2>, typename result_of::size1< T >::type >,

        // row_major <-> column_major
        std::pair<
            tag::data_order,
            typename std::conditional<
                is_column_major< T >::value, 
                tag::row_major,
                tag::column_major >::type
        >,

        // Conjugate transform (if passed by template argument)
        Conj,

        // If T has a linear array:
        // stride1 <-> stride2
        typename std::conditional< has_linear_array< T >::value,
            std::pair< tag::stride_type<1>, typename result_of::stride2< T >::type >,
            std::false_type
        >::type,
        typename std::conditional< has_linear_array< T >::value,
            std::pair< tag::stride_type<2>, typename result_of::stride1< T >::type >,
            std::false_type
        >::type,

        // If a data_side tag is present:
        // upper <-> lower
        typename std::conditional<
            property_map_has_key< prop_of_T, tag::data_side >::value,
            typename std::conditional<
                std::is_same<
                    typename property_map_at< prop_of_T, tag::data_side >::type,
                    tag::upper
                >::value,
                std::pair< tag::data_side, tag::lower >,
                std::pair< tag::data_side, tag::upper >
            >::type,
            std::false_type
        >::type

    >::type property_map;

    // Flip size1/size2
    static typename result_of::size2< T >::type size1( const Id& id ) {
        return bindings::size2( id.get() );
    }

    static typename result_of::size1< T >::type size2( const Id& id ) {
        return bindings::size1( id.get() );
    }

    // Value array access
    static typename result_of::begin_value< T >::type begin_value( Id& id ) {
        return bindings::begin_value( id.get() );
    }

    static typename result_of::end_value< T >::type end_value( Id& id ) {
        return bindings::end_value( id.get() );
    }

    // Linear array storage transpose
    // Flip stride1/stride2
    static typename result_of::stride2< T >::type stride1( const Id& id ) {
        return bindings::stride2( id.get() );
    }

    static typename result_of::stride1< T >::type stride2( const Id& id ) {
        return bindings::stride1( id.get() );
    }

};

} // namespace detail

namespace result_of {

template< typename T >
struct trans {
    typedef detail::trans_wrapper<T, std::false_type> type;
};

}

template< typename T >
typename result_of::trans<T>::type const trans( T& underlying ) {
    return detail::trans_wrapper<T, std::false_type>( underlying );
}

template< typename T >
typename result_of::trans<const T>::type const trans( const T& underlying ) {
    return detail::trans_wrapper<const T, std::false_type>( underlying );
}

} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
