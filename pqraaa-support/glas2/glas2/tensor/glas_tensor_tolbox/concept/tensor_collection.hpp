//  (C) Copyright Karl Meerbergen (2011).
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_tensor_concept_tensor_collection_hpp
#define glas_toolbox_tensor_concept_tensor_collection_hpp

#include <glas/concept/implicitly_defined.hpp>
#include <glas/concept/pure.hpp>
#include <glas/concept/collection.hpp>
#include <glas/toolbox/tensor/concept/dense_tensor_collection.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/utility/enable_if.hpp>


namespace glas { 

template <typename X, typename EnableIf=void>
struct TensorCollection
: boost::mpl::or_< glas::implicitly_defined< TensorCollection<X> >
                 , glas::DenseTensorCollection<X>
                 >
{};

template <typename X>
struct TensorCollection<X const>
: TensorCollection<X>
{};

template <typename X>
struct TensorCollection<X&>
: TensorCollection<X>
{};

} // namespace glas

// Specializations for parent concepts

namespace glas { 

template <typename X>
struct implicitly_defined<Collection<X>, typename boost::enable_if< boost::mpl::and_< glas::TensorCollection<X>
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
struct pure< glas::Collection<X>
           , typename boost::enable_if< glas::TensorCollection<X> >::type >
: boost::mpl::false_
{};

} // namespace glas

namespace glas {

template <typename X>
struct pure< glas::TensorCollection<X>
           , typename boost::enable_if< glas::DenseTensorCollection<X> >::type >
: boost::mpl::false_
{};

} // namespace glas

namespace glas { 

} // namespace glas

#endif
