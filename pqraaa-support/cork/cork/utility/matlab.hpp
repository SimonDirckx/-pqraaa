#ifndef cork_utility_matlab_hpp
#define cork_utility_matlab_hpp

#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <utility>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <complex>
#include <limits>

namespace CORK {

  template <typename T>
  typename std::enable_if< std::is_arithmetic<T>::value, std::string>::type matlab( T const& v ) {
    std::ostringstream s ;
    s << std::setprecision( std::numeric_limits<T>::digits10 ) << v ;
    return s.str()  ;
  }

  template <typename T>
  std::string matlab( std::complex<T> const& v ) {
    if (v.imag()>=0) 
      return matlab(v.real()) + "+" + matlab(v.imag()) + "i" ;
    else
      return matlab(v.real()) + "-" + matlab(-v.imag()) + "i" ;
  }

  template <typename T>
  struct matlab_vector {
    matlab_vector( T const& t, std::string const& name )
    : t_( t )
    , name( name )
    {}

    T const& t_ ;
    std::string const& name ;
  } ;

  template <typename T>
  std::ostream& operator<<( std::ostream& s, matlab_vector<T> const& t ) {
    s << t.name << " = [" ;
    if (t.t_.size()>0) s << matlab(t.t_(0)) ;
    for (typename T::size_type i=1; i<t.t_.size(); ++i) {
      s << "," << matlab(t.t_(i)) ;
    }
    s << "];\n" ;
    return s ;
  }

  template <typename T>
  struct matlab_matrix {
    matlab_matrix( T const& t, std::string const& name )
    : t_( t )
    , name( name )
    {}

    T const& t_ ;
    std::string const& name ;
  } ;

  template <typename T>
  std::ostream& operator<<( std::ostream& s, matlab_matrix<T> const& t ) {
    s << t.name << " = [" ;
    for (typename T::size_type i=0; i<t.t_.num_rows(); ++i) {
      if (t.t_.num_columns()>0) s << matlab(t.t_(i,0)) ;
      for (typename T::size_type j=1; j<t.t_.num_columns(); ++j) {
        s << "," << matlab(t.t_(i,j)) ;
      }
      s << "\n" ;
    }
    s << "];\n" ;
    return s ;
  }

  template <typename T>
  typename std::enable_if< glas2::is<glas2::DenseVector,T>::value, matlab_vector<T> >::type matlab( T const& t, std::string const& str="ans") {
    return matlab_vector<T>( t, str ) ;
  }

  template <typename T>
  typename std::enable_if< glas2::is<glas2::DenseMatrix,T>::value, matlab_matrix<T> >::type matlab( T const& t, std::string const& str="ans") {
    return matlab_matrix<T>( t, str ) ;
  }

} // namespace CORK

#endif
