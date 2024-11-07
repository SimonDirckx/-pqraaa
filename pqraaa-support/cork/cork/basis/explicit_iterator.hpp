//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the CORK Software
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_explicit_iterator_hpp
#define cork_basis_explicit_iterator_hpp

#include <glas2/vector.hpp>

namespace CORK { namespace basis {

  template <typename T>
  class explicit_iterator
  {
    public:
      typedef T   value_type ;
      typedef int size_type ;

    private:
      typedef glas2::shared_vector<value_type> temp_type ;

    public:
      template <typename Basis>
      explicit explicit_iterator( Basis const& basis, value_type const& arg )
      : temp_( basis.num_terms() )
      , index_( 0 )
      {
        basis.evaluate( arg, temp_ ) ;
      }

    public:
      value_type operator*() const { return temp_(index_) ; }

      // ++it
      explicit_iterator& operator++() {
        ++index_ ;
        return *this ;
      }

      // it++ ;
      explicit_iterator operator++(int) {
        explicit_iterator v = *this ;
        ++index_ ;
        return v ;
      }

    private:
      temp_type       temp_ ;
      size_type       index_ ;
  } ; // class explicit_iterator

} } // namespace CORK::basis

#endif
