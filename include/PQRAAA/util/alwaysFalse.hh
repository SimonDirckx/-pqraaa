#ifndef ALWAYS_FALSE_HH
#define ALWAYS_FALSE_HH
//
// Project     : Build
// File        : alwaysFalse.hh
// Description : implements constexpr returning always false, triggered only when the template is instantiated
// Author      : Kobe Bruyninckx
// Copyright   : KU Leuven Dept. CS 2023. 
//

namespace Util {

template <typename T>
constexpr bool always_false = false ;

} // namespace Util


#endif // ALWAYS_FALSE_HH