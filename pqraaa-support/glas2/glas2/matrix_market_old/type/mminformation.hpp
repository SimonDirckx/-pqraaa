//  (C) Copyright Karl Meerbergen 2009.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_matrix_market_type_mminformation_hpp
#define glas2_matrix_market_type_mminformation_hpp

#include <string>
#include <sstream>
#include <fstream>

namespace glas2 { namespace matrix_market {

  enum rep_type { ARRAY, COORDINATE } ;
  enum field_type { REAL, COMPLEX, INTEGER, PATTERN } ;
  enum symmetry_type { SYMMETRIC, HERMITIAN, SKEW_SYMMETRIC, GENERAL } ;

  class mminformation {
    public:
      mminformation()
      {}

      mminformation( std::istream& is ) {
        read( is ) ;
      }

    public:
      mminformation& operator=( mminformation const& that ) {
        num_rows_ = that.num_rows_ ;
        num_columns_ = that.num_columns_ ;
        num_nz_ = that.num_nz_ ;
        rep_ = that.rep_ ;
        field_= that.field_ ;
        symmetry_= that.symmetry_ ;
        
        return *this ;
      }

    private:
      inline void read( std::string const& file_name ) {
        std::ifstream ifs( file_name.c_str() ) ;
        read( ifs ) ;
      }

      inline std::istream& read( std::istream& is ) {
        // Read header
        std::string str ;
        std::getline( is, str, ' ' ) ;

        if ( str!="%%MatrixMarket" ) {
           throw std::string( "Not a Matrix Market file: I got on first line: ") + str ;
        }

        std::getline( is, str, ' ' ) ;
        if ( str!="matrix" && str!="mtx" ) {
           throw std::string( "Not a Matrix Market file for a matrix" ) ;
        }

        // Read matrix type
        std::getline( is, str, ' ' ) ;
        if (str=="array") rep() = ARRAY ;
        else if (str=="coordinate") rep() = COORDINATE ;
        else throw std::string( "Invalid rep" ) ;

        std::getline( is, str, ' ' ) ;
        if (str=="real") field() = REAL ;
        else if (str=="complex") field() = COMPLEX ;
        else if (str=="integer") field() = INTEGER ;
        else if (str=="pattern") field() = PATTERN ;
        else throw std::string( "Invalid value" ) ;

        std::getline( is, str ) ;
        if (str=="symmetric") symmetry() = SYMMETRIC ;
        else if (str=="skew-symmetric") symmetry() = SKEW_SYMMETRIC ;
        else if (str=="general") symmetry() = GENERAL ;
        else if (str=="hermitian") symmetry() = HERMITIAN ;
        else throw std::string( "Invalid symmetry" ) ;

        std::getline( is, str ) ;
        while ('%'==str[0]) std::getline(is, str) ;

        std::istringstream istr( str ) ;
        if (rep()==COORDINATE) istr >> num_rows() >> num_columns() >> num_nz() ;
        else {
          istr >> num_rows() >> num_columns() ;
          num_nz() = num_rows() * num_columns() ;
        }

        return is ;
      }

    public:
      typedef int size_type ;
      inline int num_rows() const { return num_rows_ ; }
      inline int& num_rows() { return num_rows_ ; }

      inline int num_columns() const { return num_columns_ ; }
      inline int& num_columns() { return num_columns_ ; }

      inline int num_nz() const { return num_nz_ ; }
      inline int& num_nz() { return num_nz_ ; }

      inline rep_type rep() const { return rep_ ; }
      inline rep_type& rep() { return rep_ ; }

      inline field_type field() const { return field_ ; }
      inline field_type& field() { return field_ ; }

      inline symmetry_type symmetry() const { return symmetry_ ; }
      inline symmetry_type& symmetry() { return symmetry_ ; }

    private:
      int           num_rows_ ;
      int           num_columns_ ;
      int           num_nz_ ;
      rep_type      rep_ ;
      field_type    field_ ;
      symmetry_type symmetry_ ;
  } ;

} } // namespace glas2::matrix_market

#endif
