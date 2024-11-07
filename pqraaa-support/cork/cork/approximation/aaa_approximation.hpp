//  (C) Copyright Karl Meerbergen 2018.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_aaa_approximation_hpp
#define cork_approximation_aaa_approximation_hpp

#include <glas2/vector.hpp>
#include <cmath>
#include <limits>

namespace CORK { namespace approximation {

  template <typename T>
  class aaa_approximation {
  public:
    typedef T                       value_type ;
    typedef decltype(std::abs(T())) real_type ;

  private:
    typedef glas2::vector<T> vector_type ;

  public:
    typedef typename vector_type::size_type size_type ;

  public:
    aaa_approximation(size_type n_buf = 100)
    : n_(0)
    , weights_(n_buf)
    , nodes_(n_buf)
    , coefficients_(n_buf)
    {}

  public:
    size_type& n() { return n_ ; }
    size_type n() const { return n_ ; }

  public:
    auto nodes() const { return nodes_( glas2::range(0,n_) ) ; }
    auto weights() const { return weights_( glas2::range(0,n_) ) ; }
    auto coefficients() const { return coefficients_( glas2::range(0,n_) ) ; }
    real_type& error() { return error_ ;}
    real_type const& error() const { return error_ ;}

  public:
    void reset( size_type n ) {
      n_ = n ;
      if (nodes_.size()!=n_) {
        nodes_.resize( n_ ) ;
        weights_.resize( n_ ) ;
        coefficients_.resize( n_ ) ;
      }
    } // reset()

  public:
    void add_node( T const& weight, T const& node, T const& coefficient ) {
      if (nodes_.size()==n_) {
        vector_type temp( n_ ) ;
        temp = nodes_(glas2::range(0,n_)) ; nodes_.resize( n_+10 ) ; nodes_(glas2::range(0,n_)) = temp ;
        temp = weights_(glas2::range(0,n_)) ; weights_.resize( n_+10 ) ; weights_(glas2::range(0,n_)) = temp ;
        temp = coefficients_(glas2::range(0,n_)) ; coefficients_.resize( n_+10 ) ; coefficients_(glas2::range(0,n_)) = temp ;
      }
      weights_(n_) = weight ;
      nodes_(n_) = node ;
      coefficients_(n_) = coefficient ;
      ++n_ ;
    } // add_node()
    
    T eval( T const& z ) const{
      T sum = 0.0 ;
      T y = 0.0 ;
      for (size_type i=0; i<n_; ++i) {
        if (z==nodes_(i)) return coefficients_(i) ;
        y += weights_(i) * coefficients_(i) / (z - nodes_(i) ) ;
        sum += weights_(i) / (z - nodes_(i) ) ;
      }

      return y / sum ;
    } // eval()
    T operator ()( T const& z )const{
      return this->eval(z);
    }

  private:
    size_type        n_ ;
    glas2::shared_vector<T> weights_ ;
    glas2::shared_vector<T> nodes_ ;
    glas2::shared_vector<T> coefficients_ ;
    real_type               error_ ;
  } ;

} } // namespace CORK::approximation

#endif
