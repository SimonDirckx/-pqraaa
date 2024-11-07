#ifndef _CONSTANTS_HH
#define _CONSTANTS_HH
//
// Project     : Constants
// File        : constants.hh
// Description : Defines constants used in the code
// Author      : Kobe Bruyninckx
// Copyright   : KU Leuven Dept. CS 2023. 
//

#include <limits>
#include <complex>

#include <BEACHpack/alias.hh>

namespace Constants {
    
using namespace Alias ;
using namespace std::complex_literals ;

#ifndef M_PIl
#define M_PIl 3.141592653589793238462643383279502884L  // define M_PIl if not defined
#endif

// General constants
    // Value of PI for different floating point precisions (std::numbers::pi only supported in C++20)
template<typename _Tp>
using _Enable_if_floating = std::enable_if_t<std::is_floating_point<_Tp>::value, _Tp>;
template < typename Tval >
inline constexpr Tval pi = _Enable_if_floating<Tval>(M_PIl) ;
    // Imaginary unit for different floating point precisions
template < typename Tval >
inline constexpr std::complex<Tval> i = std::complex<_Enable_if_floating<Tval>>(1il) ;

// LRCross constants
    // Default starting size of low-rank crossables
static const size_t default_kstart = 10 ;
    // Safety factor for recompression of low-rank crossables
template < typename Tval >
static const Tval lr_recompression_sf = 1. ; // Safety factor should be factored in when calling the recompression of LRCrossables!

// sphereSolutions constants
    // Minimum degree of the spherical basis functions in the sphere solutions
static const Tidx sph_degree_min = 10 ;
    // Maximum degree of the spherical basis functions in the sphere solutions
static const Tidx sph_degree_max = 100 ;
    // Threshold in the stopping criterion of the sphere solutions
template < typename Tval >
static const Tval sph_sol_threshold = 1e-4*std::numeric_limits<Tval>::epsilon() ;   
    // Number of consecutive times the stopping criterion should be satisfied in the sphere solutions
static const Tidx sph_sol_nb_success = 5 ;

// UHMatrix constants
    // Default safety factor for UH-matrix cluster basis computation
template < typename Tval >
static const Tval default_uhmat_sf = 0.8 ;

} // namespace Constants

#endif //_CONSTANTS_HH

