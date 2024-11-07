#ifndef glas2_bindings_detail_data_order_hpp
#define glas2_bindings_detail_data_order_hpp

#include <boost/numeric/bindings/detail/convert_to.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template<>
struct convert_to< bindings::tag::data_order, ::glas2::row_major > {
  typedef bindings::tag::row_major type;
};

template<>
struct convert_to< bindings::tag::data_order, ::glas2::column_major > {
  typedef bindings::tag::column_major type;
};


} // namespace detail
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
