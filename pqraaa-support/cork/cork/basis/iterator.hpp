//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the CORK Software
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_iterator_hpp
#define cork_basis_iterator_hpp

namespace CORK { namespace basis {

  template <typename T, typename Increment>
  class iterator
  {
    public:
      typedef T value_type ;

    public:
      explicit iterator()
      {}

      explicit iterator( Increment const& increment, value_type const& value=1.0 )
      : value_( value )
      , increment_( increment )
      {}

    public:
      value_type operator*() const { return value_ ; }

      // ++it
      iterator& operator++() {
        value_ = increment_( value_ ) ;
        return *this ;
      }

      // it++ ;
      iterator operator++(int) {
        iterator v = value_ ;
        value_ = increment_( value_ ) ;
        return v ;
      }

    private:
      value_type       value_ ;
      Increment        increment_ ;
  } ; // class iterator

} } // namespace CORK::basis

#endif
