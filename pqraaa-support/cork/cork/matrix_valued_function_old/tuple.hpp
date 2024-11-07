//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_matrix_valued_function_tuple_hpp
#define cork_matrix_valued_function_tuple_hpp

#include <cork/utility/pass_value.hpp>
#include <cork/utility/pass_lvalue_reference.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/utility/value_type_for.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace matrix_valued_function {

  template <class... NEPs>
  class tuple
  {
    public:
      template <typename ValueType>
      using value_type_for = typename std::common_type< typename coefficient_matrices_type::template value_type_for<ValueType>
                                                      , typename basis_type::template value_type_for<ValueType>
                                                      , typename coefficient_matrices_type::template value_type_for< typename basis_type::template value_type_for<ValueType> >
                                                      >::type ;

      // For testing the shift
      template <typename ValueType>
      using has_value_type_for = std::bool_constant< basis_type::template has_value_type_for<ValueType>::value
                                                   && coefficient_matrices_type::template has_value_type_for< ValueType >::value
                                                   && coefficient_matrices_type::template has_value_type_for< typename basis_type::template value_type_for<ValueType> >::value
                                                   > ;

    public:
      tuple( std::tuple<... NEPs> const& neps )
      : neps_( neps )
      {}

      static_assert( std::tuple_size< std::tuple<... NEPs> >::value>0, "CORK::matrix_valued_function::tuple: tuple is empty" ) ;

    public:
      auto size() const {
        assert( tuple<
        return std::get<0>(neps_).size() ;
      }

      auto num_terms() const {
        return basis_.num_terms() ;
      } // num_terms()

      basis_type const& basis() const { return basis_; }
      coefficient_matrices_type const& coefficient_matrices() const { return matrices_ ; }

    public:
      template <typename Shift, typename X, typename W>
      void multiply_add( Shift const& shift, X const& x, W w ) const {
        glas2::vector< typename basis_type::template value_type_for<Shift> > coefs( basis_.num_terms() ) ;
        basis_.evaluate( shift, coefs ) ;

        // Can lose efficiency. To be fixed.
/*#ifndef NDEBUG
// can be that the same multiplication is done multiple times
std::cerr << "tuple::tuple() : Can lose efficiency here" << std::endl ;
#endif*/
        for (typename decltype(coefs)::size_type i=0; i<coefs.size(); ++i)
          matrices_.multiply_add( i, coefs(i)*x, w ) ;
      } // multiply_add()

    private:
      Basis               basis_ ;
      CoefficientMatrices matrices_ ;
  } ; // class tuple


} } // namespace CORK::matrix_valued_function

#endif
