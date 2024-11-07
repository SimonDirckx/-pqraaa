//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_barycentric_rational_real_strong_hpp
#define cork_basis_barycentric_rational_real_strong_hpp

#include <cork/basis/barycentric_rational_real.hpp>

namespace CORK { namespace basis {

  template <typename Weights, typename Nodes, typename I=typename std::decay< Weights >::type::size_type>
  class barycentric_rational_real_strong
  : public barycentric_rational_real< Weights, Nodes, I >
  {
    public:
      explicit barycentric_rational_real_strong( Weights weights, Nodes nodes )
      : barycentric_rational_real< Weights, Nodes, I >( weights, nodes )
      {
        assert( nodes.size()>1 ) ;
      }

    public:
      I num_terms() const {
        return barycentric_rational_real< Weights, Nodes, I >::num_terms() ;
      }

      template <typename ShiftValueType, typename FunctionValues>
      void evaluate( ShiftValueType const& arg, FunctionValues values ) const {
        assert( values.size()==num_terms() ) ;
        barycentric_rational_real< Weights, Nodes, I >::evaluate( arg, values ) ;
      } // evaluate()
  } ; // class barycentric_rational_real

} } // namespace CORK::basis

#endif
