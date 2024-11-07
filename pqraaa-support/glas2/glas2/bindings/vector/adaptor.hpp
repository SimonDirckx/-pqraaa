#ifndef glas2_bindings_vector_adaptor_hpp
#define glas2_bindings_vector_adaptor_hpp

#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename E, typename Id >
struct adaptor< E, Id
              , typename std::enable_if< glas2::is<glas2::ContiguousDenseVector,E>::value >::type
              > {

    typedef typename copy_const< Id, typename E::value_type >::type value_type;
    /*typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::vector >,
        mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::stride_type<1>, tag::contiguous >
    > property_map;
    */
    typedef std::tuple< std::pair< tag::value_type, value_type >
                      , std::pair< tag::entity, tag::vector >
                      , std::pair< tag::size_type<1>, std::ptrdiff_t >
                      , std::pair< tag::data_structure, tag::linear_array >
                      , std::pair< tag::stride_type<1>, tag::contiguous >
                      > property_map ;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size();
    }

    static value_type* begin_value( Id& id ) {
        return id.ptr();
    }

    static value_type* end_value( Id& id ) {
        return id.ptr() + id.size();
    }

};


template< typename T, typename S, typename Id >
struct adaptor< glas2::strided_vector<T,S>, Id > {

    typedef typename copy_const< Id, T >::type value_type;
    typedef S                                  size_type1;
    typedef S                                  stride_type1;

    /*
    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::vector >,
        mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::stride_type<1>, std::ptrdiff_t >
    > property_map;
    */
    typedef std::tuple< std::pair< tag::value_type, value_type >
                      , std::pair< tag::entity, tag::vector >
                      , std::pair< tag::size_type<1>, std::ptrdiff_t >
                      , std::pair< tag::data_structure, tag::linear_array >
                      , std::pair< tag::stride_type<1>, std::ptrdiff_t >
                      > property_map ;

    // Only called in case of dynamic sizes
    static std::ptrdiff_t size1( const Id& id ) {
        return id.size();
    }

    static value_type* begin_value( Id& id ) {
        return id.ptr();
    }

    static value_type* end_value( Id& id ) {
        return id.ptr() + id.stride()*id.size() ;
    }

    // Only called in case of dynamic strides
    static std::ptrdiff_t stride1( const Id& id ) {
        return id.stride();
    }

};


} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
