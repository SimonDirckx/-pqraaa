//  (C) Copyright Karl Meerbergen 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_accumulated_partial_fractions_real_hpp
#define cork_basis4cork_accumulated_partial_fractions_real_hpp

#include <cork/utility/is_infinite.hpp>
#include <cork/basis/accumulated_partial_fractions_real.hpp>
#include <cork/basis4cork/accumulated_partial_fractions_real_matrices.hpp>
#include <cork/basis4cork/explicit_matrices.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <cork/basis4cork/evaluate.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/matrix.hpp>
#include <glas2/bindings/vector.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getrs.hpp>
#include <cassert>
#include <string>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename Weights, typename Nodes, typename I>
  class basis4CORK< basis::accumulated_partial_fractions_real<Weights,Nodes,I> >
  : public explicit_matrices< decltype( std::abs(typename basis::accumulated_partial_fractions_real<Weights,Nodes,I>::data_value_type()) ) >
  {
    public:
      typedef typename std::decay<Weights>::type weights_type ;
      typedef typename std::decay<Nodes>::type   poles_type ;
      typedef I                                  size_type ;

    public:
      explicit basis4CORK( basis::accumulated_partial_fractions_real<Weights,Nodes,I> const& basis )
      : explicit_matrices< decltype( std::abs(typename basis::accumulated_partial_fractions_real<Weights,Nodes,I>::data_value_type()) ) >( basis.num_terms() )
      , weights_( basis.weights() )
      , poles_( basis.poles() )
      {
        accumulated_partial_fractions_real_matrices( weights_, poles_, basis.num_real(), this->M(), this->N() ) ;
      }

    public:
      // Basis4CORK
      int num_terms() const { return poles_.size()+1 ; }

      typename std::add_lvalue_reference< typename std::add_const<Weights>::type >::type weights() const& {
        return weights_ ;
      }

      typename std::add_lvalue_reference< typename std::add_const<Nodes>::type >::type poles() const& {
        return poles_ ;
      }

    public:
      template <typename Coefs>
      void evaluate( Coefs coefs ) const {
        CORK::basis4CORK::evaluate( *this, coefs ) ;
      }

    private:
      weights_type const& weights_ ;
      poles_type const&   poles_ ;
  } ; // accumulated_partial_fractions_real

} } // namespace CORK::basis4cork

#endif
