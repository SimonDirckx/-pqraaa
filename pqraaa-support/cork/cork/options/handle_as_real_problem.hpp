//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef CORK_OPTIONS_HANDLE_AS_REAL_PROBLEM__HPP
#define CORK_OPTIONS_HANDLE_AS_REAL_PROBLEM__HPP

#include<tuple>
#include<type_traits>

namespace CORK { namespace options {

  class handle_as_real_problem_tag
  {
    public:
      inline handle_as_real_problem_tag()
      {}

      template <typename ...Ts>
      inline handle_as_real_problem_tag( std::tuple<Ts...> const& options)
      {}

      inline std::false_type value() const { return std::false_type() ; }
  } ;

  class is_real_problem
  : public handle_as_real_problem_tag
  {
    public:
      is_real_problem()
      {}

      template <typename ...Ts>
      is_real_problem( std::tuple<Ts...> const& options)
      {}

      inline std::true_type value() const { return std::true_type() ; }
  } ;

  inline auto handle_as_real_problem() { return is_real_problem() ; }

} } // CORK::options

#endif
