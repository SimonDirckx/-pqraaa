//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_sv_aaa_approximation_hpp
#define cork_approximation_sv_aaa_approximation_hpp

#include <cork/approximation/info.hpp>
#include <cork/approximation/sv_aaa_triple.hpp>
#include <cmath>
#include <limits>
#include <cstring>
#include <vector>

namespace CORK { namespace approximation {

  template <typename Arg, typename Result, typename TF=Result>
  class SV_AAA_approximation
  : private SV_AAA_triple<Arg,Result,TF> {
  public:
    typedef typename SV_AAA_triple<Arg,Result,TF>::value_type value_type ;
    typedef decltype(std::abs(value_type()))      real_type ;

  public:
    typedef int size_type ;

  public:
    SV_AAA_approximation(size_type n_fun, size_type n_buf = 100, std::string const& name="SV-AAA")
    : SV_AAA_triple<Arg,Result,TF>( n_fun, n_buf )
    , info_( name )
    {}

  public:
    auto nodes() const { return this->nodes_( glas2::range(0,this->n()) ) ; }
    auto weights() const { return this->weights_( glas2::range(0,this->n()) ) ; }
    auto coefficients() const { return this->coefficients_( glas2::range(0,this->n()), glas2::all() ) ; }

  public:
    void reset( size_type n ) {
      assert( n<=this->n() ) ;
      this->n() = n ;
      static_cast< SV_AAA_triple<Arg,Result,TF>& >(*this).reset( n ) ;
    } // reset()

    int n() const {return info_.size ; }
    int& n() {return info_.size ; }

  public:
    template <typename Coefficient>
    void add_node( Result const& weight, Arg const& node, Coefficient const& coefficient ) {
      if (this->nodes_.size()==info_.size) {
        {
          glas2::vector<Arg>  temp( info_.size ) ;
          temp = this->nodes_(glas2::range(0,info_.size)) ; this->nodes_.resize( info_.size+10 ) ; this->nodes_(glas2::range(0,info_.size)) = temp ;
        }
        {
          glas2::vector<Result>  temp( info_.size ) ;
          temp = this->weights_(glas2::range(0,info_.size)) ; this->weights_.resize( info_.size+10 ) ; this->weights_(glas2::range(0,info_.size)) = temp ;
        }
        {
          glas2::matrix<TF> temp( info_.size, this->coefficients_.num_columns() ) ;
          temp = this->coefficients_(glas2::range(0,info_.size), glas2::all()) ; this->coefficients_.resize( info_.size+10, this->coefficients_.num_columns() ) ; this->coefficients_(glas2::range(0,info_.size), glas2::all()) = temp ;
        }
      }
      this->weights_(info_.size) = weight ;
      this->nodes_(info_.size) = node ;
      this->coefficients_(info_.size, glas2::all() ) = coefficient ;
      ++info_.size ;
    } // add_node()

    template <typename Vector>
    void eval( Arg const& z, Vector result ) {
      Result sum = 0.0 ;
      fill( result, 0.0 ) ;
      for (size_type i=0; i<info_.size; ++i) {
        if (z==this->nodes_(i)) {
          result = this->coefficients_(i, glas2::all()) ;
          return ;
        }
        if (std::abs(this->nodes_(i))!=std::numeric_limits<real_type>::infinity()) {
          result += this->weights_(i) * this->coefficients_(i, glas2::all()) / (z - this->nodes_(i) ) ;
          sum += this->weights_(i) / (z - this->nodes_(i) ) ;
        } else {
          result += this->weights_(i) * this->coefficients_(i, glas2::all()) ;
          sum += this->weights_(i) ;
        }
      }
      if (info_.size>0) result /= sum ;
    } // eval()

  private:
    approximation::info<double> info_ ; // Can use real_type instead

  public:
    approximation::info<real_type> const& info() const { return info_ ; }
    std::vector< std::string >& warnings() { return info_.warnings ; } ;
    real_type& error() { return info_.error ; } ;
    real_type const& error() const { return info_.error ; } ;
  } ; // class SV_AAA_approximation

} } // namespace CORK::approximation

#endif
