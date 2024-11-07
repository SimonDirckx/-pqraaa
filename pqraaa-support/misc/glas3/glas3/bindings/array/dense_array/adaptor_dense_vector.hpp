//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_bindings_array_dense_array_adaptor_dense_vector_hpp
#define glas3_bindings_array_dense_array_adaptor_dense_vector_hpp

#include <glas3/array/dense_array/concept/contiguous_dense_array.hpp>
#include <glas3/concept/is.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename E, typename Id >
struct adaptor< E, Id, typename std::enable_if< glas3::is<glas3::ContiguousDenseVector, E>::value >::type > {

    typedef typename copy_const< Id, typename E::value_type >::type value_type;
    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::vector >,
        mpl::pair< tag::size_type<1>, typename E::size_type >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::stride_type<1>, tag::contiguous >
    > property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size() ;
    }

    static value_type* begin_value( Id& id ) {
        return id.data_ptr().get() ;
    }

    static value_type* end_value( Id& id ) {
        return id.data_ptr().get() + id.size() ;
    }

};


//template< typename T, typename S, typename Id >
//struct adaptor< glas2::strided_vector<T,S>, Id > {
//
//    typedef typename copy_const< Id, T >::type value_type;
//    typedef S                                  size_type1;
//    typedef S                                  stride_type1;
//
//    typedef mpl::map<
//        mpl::pair< tag::value_type, value_type >,
//        mpl::pair< tag::entity, tag::vector >,
//        mpl::pair< tag::size_type<1>, size_type1 >,
//        mpl::pair< tag::data_structure, tag::linear_array >,
//        mpl::pair< tag::stride_type<1>, stride_type1 >
//    > property_map;
//
//    // Only called in case of dynamic sizes
//    static std::ptrdiff_t size1( const Id& id ) {
//        return id.size();
//    }
//
//    static value_type* begin_value( Id& id ) {
//        return id.ptr();
//    }
//
//    static value_type* end_value( Id& id ) {
//        return id.ptr() + id.stride()*id.size() ;
//    }
//
//    // Only called in case of dynamic strides
//    static std::ptrdiff_t stride1( const Id& id ) {
//        return id.stride();
//    }
//
//};


} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
