//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_barycentric_rational_strong_hpp
#define cork_basis_barycentric_rational_strong_hpp

#include <cork/basis/barycentric_rational.hpp>

namespace CORK { namespace basis {

  template <typename Weights, typename Nodes, typename I=typename std::decay< Weights >::type::size_type>
  class barycentric_rational_strong
  : public barycentric_rational< Weights, Nodes, I>
  {
    private:
      typedef barycentric_rational< Weights, Nodes, I> base_type ;

    public:
      explicit barycentric_rational_strong( Weights weights, Nodes nodes )
      : base_type( weights, nodes )
      {
        assert( weights.size()>1 ) ;
        assert( nodes.size()==weights.size() ) ;
      } // barycentric_rational

    public:
      typedef typename base_type::minimal_value_type minimal_value_type ;

/*    public:
      template <typename ValueType>
      using value_type_for = typename base_type::template value_type_for<ValueType> ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< T, inner_value_type > > ;
*/
    public:
      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        static_cast<base_type const&>( *this ).evaluate( arg, values ) ;
      } // evaluate()
  } ; // class barycentric_rational


} } // namespace CORK::basis

#endif
