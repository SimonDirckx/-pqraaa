//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_pf_approximation_hpp
#define cork_approximation_pf_approximation_hpp

#include <glas2/vector.hpp>
#include <cmath>
#include <limits>

namespace CORK { namespace approximation {

  template <typename T>
  class pf_approximation {
  public:
    typedef T                       value_type ;
    typedef decltype(std::abs(T())) real_type ;

  public:
    typedef glas2::shared_vector<T> vector_type ;
    typedef glas2::shared_matrix<T> matrix_type ;

  public:
    typedef typename vector_type::size_type size_type ;

  public:
    template <typename W, typename N, typename C>
    pf_approximation( W const& weights, N const& nodes, C const& coefficients )
    : weights_(glas2::copy(weights))
    , nodes_(glas2::copy(nodes))
    , coefficients_(coefficients.num_rows(), coefficients.num_columns())
    {
      assert( weights.size()==nodes.size() ) ;
      assert( weights.size()+1==coefficients.num_rows() ) ;
      coefficients_ = coefficients ;
    }

  public:
    size_type n() const { return nodes_.size() ; }

  public:
    auto const& nodes() const { return nodes_ ; }
    auto const& weights() const { return weights_ ; }
    auto const& coefficients() const { return coefficients_ ; }

  public:
    template <typename Arg, typename Vector>
    void eval( Arg const& z, Vector result ) {
      assert( result.size() == coefficients_.num_columns() ) ;
      result = coefficients_( 0, glas2::all() ) ;
      for (size_type i=0; i<n(); ++i) {
        result += (weights_(i) / (z - nodes_(i))) * coefficients_(i+1,glas2::all()) ;
      }
    } // eval()

  private:
    vector_type weights_ ;
    vector_type nodes_ ;
    matrix_type coefficients_ ;
  } ;

} } // namespace CORK::approximation

#endif
