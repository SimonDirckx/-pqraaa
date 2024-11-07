//  (C) Copyright Karl Meerbergen 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_type_mm_reader_hpp
#define glas2_matrix_market_type_mm_reader_hpp

#include <glas2/sparse/concept/forward_coordinate_sparse_matrix.hpp>
#include <glas2/concept/concept.hpp>
#include <string>
#include <sstream>
#include <fstream>
#ifdef GLAS_COMPLEX
#include <complex>
#endif

namespace glas2 { namespace matrix_market {

  enum rep_type { ARRAY, COORDINATE } ;
  enum field_type { REAL, COMPLEX, INTEGER, PATTERN } ;
  enum symmetry_type { SYMMETRIC, HERMITIAN, SKEW_SYMMETRIC, GENERAL } ;

  template <typename T>
  struct compose {
    template <typename T2>
    static T apply ( T2 real, T2 imag ) {
      assert( imag==0.0 ) ;
      return real ;
    }
  } ;

#ifdef GLAS_COMPLEX
  template <typename T>
  struct compose< std::complex<T> > {
    template <typename T2>
    static std::complex<T> apply ( T2 real, T2 imag ) {
      return std::complex<T>( real, imag ) ;
    }
  } ;
#endif

  template <typename T>
  class mm_reader {
    public:
      typedef T   value_type ;
      typedef int size_type ;

    public:
      mm_reader( std::istream& is )
      : s_ (is )
      {
        read() ;
      }

      operator std::istream& () { return s_ ; }

    public:
      struct iterator {
        iterator( std::istream& s, size_type index, field_type f, symmetry_type symm )
        : s_( s )
        , index_( index )
        , field_( f )
        , symm_( symm )
        , symm_switch_( false )
        {
          if (index_==0) {
            if (!symm_switch_) {
              read() ;
              if (symm_==SYMMETRIC || symm_==SKEW_SYMMETRIC) {
                if (column_!=row_) symm_switch_ = true ;
              }
            } else {
              size_type t = row_ ; row_ = column_ ; column_ = t ;
              if (symm_==SKEW_SYMMETRIC) value_ = -value_ ;
              symm_switch_ = false ;
            }
          }
        }

        void read() {
          std::getline( s_, str_ ) ;
          while ('%'==str_[0]) std::getline(s_, str_) ;
          std::stringstream s_stream( str_ ) ;
          s_stream >> row_ >> column_ ; --row_ ; --column_ ;
          if (field_==REAL) {
            decltype( std::abs(value_type()) ) value ;
            s_stream >> value ;
            value_ = value ;
          } else if (field_==COMPLEX) {
            decltype( std::abs(value_type()) ) real, imag ;
            s_stream >> real >> imag ;
            value_ = compose<value_type>::apply( real, imag ) ;
          }
        }

        void operator++() {
          if (!symm_switch_) {
            read() ;
            ++index_ ;
            if (symm_==SYMMETRIC || symm_==SKEW_SYMMETRIC) {
              if (column_!=row_) symm_switch_ = true ;
            }
          } else {
            size_type t = row_ ; row_ = column_ ; column_ = t ;
            if (symm_==SKEW_SYMMETRIC) value_ = -value_ ;
            symm_switch_ = false ;
          }
        }

        bool operator==( iterator const& that ) const {return index_==that.index_ ;}
        bool operator!=( iterator const& that ) const {return index_!=that.index_ ;}

        size_type row() const { return row_ ; }
        size_type column() const { return column_ ; }
        value_type value() const { /*std::cout << row_ << " " << column_ << " " << value_ << std::endl ;*/ return value_ ; }
        
        std::istream& s_ ;
        size_type     index_ ;
        field_type    field_ ;
        symmetry_type symm_ ;
        bool          symm_switch_ ;
        size_type     row_ ;
        size_type     column_ ;
        value_type    value_ ;
        std::string   str_ ;
      } ;

      iterator begin() const { return iterator( s_, 0, field(), symmetry() ) ; }
      iterator end() const { return iterator( s_, num_nz_, field(), symmetry() ) ; }

    private:
      inline void read() {
        // Read header
        std::string str ;
        std::getline( s_, str, ' ' ) ;

        if ( str!="%%MatrixMarket" ) {
           throw std::string( "Not a Matrix Market file: I got on first line: ") + str ;
        }

        std::getline( s_, str, ' ' ) ;
        if ( str!="matrix" && str!="mtx" ) {
           throw std::string( "Not a Matrix Market file for a matrix" ) ;
        }

        // Read matrix type
        std::getline( s_, str, ' ' ) ;
        if (str=="array") rep_ = ARRAY ;
        else if (str=="coordinate") rep_ = COORDINATE ;
        else throw std::string( "Invalid rep" ) ;

        std::getline( s_, str, ' ' ) ;
        if (str=="real") field_ = REAL ;
        else if (str=="complex") field_ = COMPLEX ;
        else if (str=="integer") field_ = INTEGER ;
        else if (str=="pattern") field_ = PATTERN ;
        else throw std::string( "Invalid value" ) ;

        std::getline( s_, str ) ;
        if (str=="symmetric") symmetry_ = SYMMETRIC ;
        else if (str=="skew-symmetric") symmetry_ = SKEW_SYMMETRIC ;
        else if (str=="general") symmetry_ = GENERAL ;
//        else if (str=="hermitian") symmetry_ = HERMITIAN ;
        else throw std::string( "Invalid symmetry for GLAS2" ) ;

        std::getline( s_, str ) ;
        while ('%'==str[0]) std::getline(s_, str) ;

        std::istringstream istr( str ) ;
        if (rep()==COORDINATE) istr >> num_rows_ >> num_columns_ >> num_nz_ ;
        else {
          istr >> num_rows_ >> num_columns_ ;
          num_nz_ = num_rows_ * num_columns_ ;
        }
      }

    public:
      inline size_type num_rows() const { return num_rows_ ; }
      inline size_type num_columns() const { return num_columns_ ; }
      inline size_type num_nz() const { return num_nz_ ; }

      inline rep_type rep() const { return rep_ ; }
      inline field_type field() const { return field_ ; }
      inline symmetry_type symmetry() const { return symmetry_ ; }

      inline std::istream& stream() { return s_ ; }

    private:
      std::istream& s_ ;
      bool                  read_ ;
      int                   num_rows_ ;
      int                   num_columns_ ;
      int                   num_nz_ ;
      rep_type              rep_ ;
      field_type            field_ ;
      symmetry_type         symmetry_ ;
  } ;

} } // namespace glas2::matrix_market

namespace glas2 {

  template <typename T>
  struct concept< matrix_market::mm_reader<T> >
  : ForwardCoordinateSparseMatrix
  {} ;

} // namespace glas2

#endif
