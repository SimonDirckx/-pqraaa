//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_triple_hpp
#define cork_approximation_sv_aaa_triple_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>

namespace CORK { namespace approximation {

  template <typename Arg, typename Result, typename TF=Result>
  class SV_AAA_triple {
  public:
    typedef Result value_type ;
    typedef int    size_type ;

  public:
    SV_AAA_triple(size_type n_fun, size_type n)
    : weights_(n)
    , nodes_(n)
    , coefficients_(n, n_fun)
    {}

  public:
    size_type size() const { return nodes_.size();}

    void reset( size_type n ) {
      nodes_.resize( n ) ;
      weights_.resize( n ) ;
      coefficients_.resize( n, coefficients_.num_columns() ) ;
    } // reset()

  public:
    glas2::shared_vector<Result> weights_ ;
    glas2::shared_vector<Arg>    nodes_ ;
    glas2::shared_matrix<TF>     coefficients_ ;
  } ; // class SV_AAA_triple

} } // namespace CORK::approximation

#endif
