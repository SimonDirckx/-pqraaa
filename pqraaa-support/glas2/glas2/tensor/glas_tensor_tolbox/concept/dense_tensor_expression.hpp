//  (C) Copyright Karl Meerbergen (2011).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_concept_dense_tensor_expression_hpp
#define glas_toolbox_tensor_concept_dense_tensor_expression_hpp

#include <glas/concept/implicitly_defined.hpp>
#include <glas/concept/pure.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/utility/enable_if.hpp>


namespace glas { 

template <typename X, typename EnableIf=void>
struct DenseTensorExpression
: glas::implicitly_defined< DenseTensorExpression<X> >
{};

template <typename X>
struct DenseTensorExpression<X const>
: DenseTensorExpression<X>
{};

template <typename X>
struct DenseTensorExpression<X&>
: DenseTensorExpression<X>
{};

} // namespace glas

namespace glas { 

} // namespace glas

#endif
