//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_UBLAS_MATRIX_EXPRESSION_HPP
#define BOOST_NUMERIC_BINDINGS_UBLAS_MATRIX_EXPRESSION_HPP

#include <boost/numeric/bindings/bandwidth.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/property_map.hpp>
#include <boost/numeric/bindings/end.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <boost/mpl/replace.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename Id, typename Enable >
struct adaptor< boost::numeric::ublas::matrix_reference< T >, Id, Enable > {

    typedef typename copy_const< Id, T >::type adapted_type;
    typedef typename property_map_of< adapted_type >::type property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size1();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.size2();
    }

    static typename result_of::begin_value< adapted_type >::type begin_value( Id& id ) {
        return bindings::begin_value( id.expression() );
    }

    static typename result_of::end_value< adapted_type >::type end_value( Id& id ) {
        return bindings::end_value( id.expression() );
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return bindings::stride1( id.expression() );
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return bindings::stride2( id.expression() );
    }

    static std::ptrdiff_t bandwidth1( const Id& id ) {
        return bindings::bandwidth1( id.expression() );
    }

    static std::ptrdiff_t bandwidth2( const Id& id ) {
        return bindings::bandwidth2( id.expression() );
    }

};

template< typename T, typename U, typename Id, typename Enable >
struct adaptor< boost::numeric::ublas::matrix_unary2< T, U >, Id, Enable > {

    typedef typename copy_const< Id, T >::type adapted_type;
    typedef typename property_map_of< adapted_type >::type map;

    typedef std::tuple<
        std::pair<tag::value_type, typename property_map_at<map, tag::value_type>::type>,
        std::pair<tag::entity, typename property_map_at<map, tag::entity>::type>,
        std::pair<tag::size_type<1>, typename property_map_at<map, tag::size_type<1> >::type>,
	std::pair<tag::size_type<2>, typename property_map_at<map, tag::size_type<2> >::type>,
	std::pair<tag::data_structure, typename property_map_at<map, tag::data_structure>::type>,

        std::pair<tag::data_order,
		  typename std::conditional<
		      std::is_same<
			  typename property_map_at<map, tag::data_order>::type,
			  tag::row_major>::value,
		      tag::column_major,
		      tag::row_major
		      >::type>,

        std::pair<tag::stride_type<1>, typename property_map_at<map, tag::stride_type<2> >:: type>,
        std::pair<tag::stride_type<2>, typename property_map_at<map, tag::stride_type<1> >:: type>
	> property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size1();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.size2();
    }

    static typename result_of::begin_value< adapted_type >::type
    begin_value( Id& id ) {
        return bindings::begin_value( id.expression() );
    }

    static typename result_of::end_value< adapted_type >::type
    end_value( Id& id ) {
        return bindings::end_value( id.expression() );
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return bindings::stride2( id.expression() );
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return bindings::stride1( id.expression() );
    }

};

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
