//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_cursor_hpp
#define cork_matrix_iterator_cursor_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK {

  template <typename DegreeType, typename Schedule>
  class cursor {
    public:
      typedef DegreeType degree_type ;
      typedef Schedule   schedule_type ;

    public:
      cursor( degree_type index, schedule_type const& schedule )
          : index_( index )
          , schedule_( schedule )
          {}

    public:
      bool operator==( cursor const& that ) const { return index_ == that.index_ ; }
      bool operator!=( cursor const& that ) const { return index_ != that.index_ ; }
      bool operator<( cursor const& that ) const { return index_ < that.index_ ; }
      bool operator>( cursor const& that ) const { return index_ > that.index_ ; }
      bool operator<=( cursor const& that ) const { return index_ <= that.index_ ; }
      bool operator>=( cursor const& that ) const { return index_ >= that.index_ ; }

    public:
      cursor& operator+=( degree_type i ) {
        index_ += i ;
        return *this ;
      }

      cursor& operator-=( degree_type i ) {
        index_ -= i ;
        return *this ;
      }

      cursor& operator++() {
        ++index_ ;
        return *this ;
      }

      cursor& operator--() {
        --index_ ;
        return *this ;
      }

    public:
      degree_type index() const { return index_ ; }

      template <typename Z>
      void schedule( Z const& z ) {
        schedule_( index_, z ) ;
      }

    private:
      degree_type   index_ ;
      schedule_type schedule_ ;
  } ; // class cursor

  template <typename DegreeType, typename Schedule>
  cursor< DegreeType, Schedule > make_cursor( DegreeType index, Schedule const& schedule ) {
    return cursor< DegreeType, Schedule >( index, schedule ) ;
  }

} // namespace CORK

#endif
