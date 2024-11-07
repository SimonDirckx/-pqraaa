#ifndef cork_utility_value_type_for_hpp
#define cork_utility_value_type_for_hpp

namespace CORK {

  template <typename T, typename Class>
  struct value_type_for {
    typedef typename Class::template value_type<T> type ;
  } ;

  template <typename T, typename Class>
  struct value_type_for<T, Class&>
  : value_type_for< T, Class >
  {} ;

  template <typename T, typename Class>
  struct value_type_for<T, Class const>
  : value_type_for< T, Class >
  {} ;

} // namespace CORK

#endif
