//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_matrix_iterator_merged_cursor_hpp
#define cork_matrix_iterator_merged_cursor_hpp

#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <type_traits>

namespace CORK {

  template <typename Cursor_1, typename Cursor_2>
  class merged_cursor {
    public:
      typedef typename std::decay< Cursor_1 >::type                                                                      cursor_1_type ;
      typedef typename std::decay< Cursor_2 >::type                                                                      cursor_2_type ;
      typedef typename std::common_type< typename cursor_1_type::degree_type, typename cursor_2_type::degree_type>::type degree_type ;

    public:
      merged_cursor( Cursor_1 cursor_1, Cursor_2 cursor_2 )
      : cursor_1_( cursor_1 )
      , cursor_2_( cursor_2 )
      {}

    public:
      bool operator==( merged_cursor const& that ) const { return index_ == that.index_ ; }
      bool operator!=( merged_cursor const& that ) const { return index_ != that.index_ ; }
      bool operator<( merged_cursor const& that ) const { return index_ < that.index_ ; }
      bool operator>( merged_cursor const& that ) const { return index_ > that.index_ ; }
      bool operator<=( merged_cursor const& that ) const { return index_ <= that.index_ ; }
      bool operator>=( merged_cursor const& that ) const { return index_ >= that.index_ ; }

    public:
      merged_cursor& operator++() {
        ++index_ ;
        return *this ;
      }

      merged_cursor& operator--() {
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
      Cursor_1 cursor_1_ ;
      Cursor_2 cursor_2_ ;
      degree_type     border_ ;
  } ; // class merged_cursor

  template <typename DegreeType, typename Schedule>
  merged_cursor< DegreeType, Schedule > make_merged_cursor( DegreeType index, Schedule const& schedule ) {
    return merged_cursor< DegreeType, Schedule >( index, schedule ) ;
  }

} // namespace CORK

#endif
