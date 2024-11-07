//  (C) Copyright Karl Meerbergen (2011).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_concept_tensor_expression_hpp
#define glas_toolbox_tensor_concept_tensor_expression_hpp

#include <glas/concept/implicitly_defined.hpp>
#include <glas/concept/pure.hpp>
#include <glas/concept/expression.hpp>
#include <glas/toolbox/tensor/concept/dense_tensor_expression.hpp>
#include <glas/toolbox/tensor/concept/tensor_collection.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/utility/enable_if.hpp>


namespace glas { 

template <typename X, typename EnableIf=void>
struct TensorExpression
: boost::mpl::or_< glas::implicitly_defined< TensorExpression<X> >
                 , glas::DenseTensorExpression<X>
                 , glas::TensorCollection<X>
                 >
{};

template <typename X>
struct TensorExpression<X const>
: TensorExpression<X>
{};

template <typename X>
struct TensorExpression<X&>
: TensorExpression<X>
{};

} // namespace glas

// Specializations for parent concepts

namespace glas { 

template <typename X>
struct implicitly_defined<glas::Expression<X>, typename boost::enable_if< boost::mpl::and_< glas::TensorExpression<X>
                                                               , boost::mpl::not_< boost::is_const< X > >
                                                               , boost::mpl::not_< boost::is_reference< X > >
                                                               >
                                             >::type
                >
: boost::mpl::true_
{};

} // namespace glas

namespace glas {

template <typename X>
struct pure< ::glas::Expression<X>
           , typename boost::enable_if< glas::TensorExpression<X> >::type >
: boost::mpl::false_
{};

} // namespace glas

namespace glas {

template <typename X>
struct pure< glas::TensorExpression<X>
           , typename boost::enable_if< glas::DenseTensorExpression<X> >::type >
: boost::mpl::false_
{};

} // namespace glas

namespace glas { 

} // namespace glas

#endif
