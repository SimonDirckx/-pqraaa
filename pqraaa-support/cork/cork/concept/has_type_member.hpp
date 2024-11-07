#ifndef CORK_CONCEPT_HAS_TYPE_MEMBER_HPP
#define CORK_CONCEPT_HAS_TYPE_MEMBER_HPP

#include <type_traits>

namespace CORK {

template< class, class = void >
struct has_type_member : std::false_type { };

template< class T >
struct has_type_member<T, std::void_t<typename T::type>> : std::true_type { };

} // namespace CORK

#endif
