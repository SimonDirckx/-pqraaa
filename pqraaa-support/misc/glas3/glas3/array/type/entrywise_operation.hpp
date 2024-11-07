//	  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_type_entrywise_operation_hpp
#define glas3_array_type_entrywise_operation_hpp

namespace glas3 {

template <typename Op, typename X, typename Y, typename EnableIf=void>
class binary_operation
{} ;

template <typename Op, typename X, typename EnableIf=void>
class unary_operation
{} ;

} // namespace glas3

#endif
