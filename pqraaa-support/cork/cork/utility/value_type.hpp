#ifndef cork_utility_value_type_hpp
#define cork_utility_value_type_hpp

namespace CORK {

  template <typename Class>
  struct value_type {
    typedef typename Class::value_type type ;
  } ;

  template <typename Class>
  struct value_type<Class&>
  : value_type<Class>
  {} ;

  template <typename Class>
  struct value_type<Class const>
  : value_type<Class>
  {} ;

} // namespace CORK

#endif
