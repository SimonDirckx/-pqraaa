#ifndef glas2_bindings_vector_std_initializer_list_hpp
#define glas2_bindings_vector_std_initializer_list_hpp

#include <glas2/vector/concept/contiguous_dense_vector.hpp>
#include <glas2/concept/is.hpp>
#include <glas2/vector/type/strided_vector.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <initializer_list>
#include <type_traits>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< typename Id, typename T >
struct adaptor< std::initializer_list<T>, Id > {

    typedef typename copy_const< Id, typename std::initializer_list<T>::value_type >::type value_type;
    /*
    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::entity, tag::vector >,
        mpl::pair< tag::size_type<1>, typename std::initializer_list<T>::size_type >,
        mpl::pair< tag::data_structure, tag::linear_array >,
        mpl::pair< tag::stride_type<1>, tag::contiguous >
    > property_map;
    */
    typedef std::tuple< std::pair< tag::value_type, value_type >
                      , std::pair< tag::entity, tag::vector >
                      , std::pair< tag::size_type<1>, typename std::initializer_list<T>::size_type >
                      , std::pair< tag::data_structure, tag::linear_array >
                      , std::pair< tag::stride_type<1>, tag::contiguous >
                      > property_map ;

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
