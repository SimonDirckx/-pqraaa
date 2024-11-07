//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_approximation_vector_of_functions_hpp
#define cork_approximation_vector_of_functions_hpp

#include <cassert>
#include <type_traits>
#include <vector>
#include <functional>

namespace CORK { namespace approximation {

  template <typename T>
  class vector_of_functions
  {
    public:
      typedef T                               value_type ;
      typedef std::function< T(T) >           function_type ;
      typedef std::vector< function_type >    vector_type ;
      typedef typename vector_type::size_type size_type ;

    public:
      explicit vector_of_functions()
      {}

      explicit vector_of_functions( vector_type const& vector )
      : vector_( vector )
      {}

    public:
      size_type num_terms() const {
        return vector_.size() ;
      } // num_terms()

      void push_back( function_type const& f ) { vector_.push_back(f) ; }
      decltype (auto) operator[]( size_type i ) const { return vector_[i] ; }

    public:
      template <typename ValueType, typename FunctionValues>
      void evaluate( ValueType const& arg, FunctionValues values ) const {
        assert( values.size() == num_terms() ) ;
        for (size_type i=0; i<num_terms(); ++i) values(i) = vector_[i]( arg ) ;
      } // evaluate

    private:
      vector_type vector_ ;

  } ; // class vector_of_functions


} } // namespace CORK::approximation

#endif
