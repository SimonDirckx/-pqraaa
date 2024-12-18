//
// Copyright (c) 2009 Rutger ter Borg
// Copyright (c) 2011 Andrey Asadchev
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_UBLAS_STORAGE_HPP
#define BOOST_NUMERIC_BINDINGS_UBLAS_STORAGE_HPP

#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename T, typename Id >
struct adaptor< T, Id, typename boost::enable_if< std::is_base_of< ublas::storage_array<T>, T> >::type > {

    typedef typename copy_const< Id, typename T::value_type >::type value_type;
    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::entity, tag::vector >,
        std::pair< tag::size_type<1>, std::ptrdiff_t >,
        std::pair< tag::data_structure, tag::linear_array >,
        std::pair< tag::stride_type<1>, tag::contiguous >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size();
    }

    static value_type* begin_value( Id& id ) {
        return id.begin();
    }

    static value_type* end_value( Id& id ) {
        return id.end();
    }

};

} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
