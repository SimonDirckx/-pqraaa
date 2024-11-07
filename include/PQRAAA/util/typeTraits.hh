#ifndef _TYPE_TRAITS_HH
#define _TYPE_TRAITS_HH
//
// Project     : Cross
// File        : crossBase.hh
// Description : implements additional type traits
// Author      : Kobe Bruyninckx
// Copyright   : KU Leuven Dept. CS 2024. 
//
#include <type_traits>

#include <BEACHpack/alias.hh>
#include <BEACHpack/bilinear/bilinearForm.hh>
#include <BEACHpack/fnSpace.hh>
#include <BEACHpack/kernel/slp.hh>
#include <BEACHpack/lowrank/compactRep.hh>
#include <BEACHpack/lowrank/LRMatrix.hh>

namespace Util {

using namespace Alias ;

/**
 * @brief Determines whether a class is crossable, providing seven member functions:
 * - getRow(i,v) fills v with the i-th row
 * - getCol(j,w) fills w with the j-th column
 * - operator()(i,j) returns the element at position (i,j)
 * - rows() returns the number of rows
 * - cols() returns the number of columns
 * - matrix() returns the full matrix, convertible to Eigen's matrix
 * - getMatrix(m) fills m with the full matrix
 * 
 * The idea of crossables is important for ACA to allow efficient computation
 * of rows and columns. This cannot be captured in Eigen's expression templates
 * as they rely on element-wise sampling of the expression.
 * 
 * @tparam TClass Type of the class
 * @tparam Tval Value type of the elements
 */
template< typename TClass, typename Tval >
struct is_crossable {
    static const bool value =   std::is_invocable<decltype(&TClass::getRow),TClass,Tidx,Vec<Tval>&>::value 
                            &&  std::is_invocable<decltype(&TClass::getCol),TClass,Tidx,Vec<Tval>&>::value 
                            &&  std::is_convertible<decltype(std::declval<TClass>().operator()(static_cast<Tidx>(0),static_cast<Tidx>(0))),Tval>::value
                            &&  std::is_convertible<decltype(std::declval<TClass>().rows()),Tidx>::value
                            &&  std::is_convertible<decltype(std::declval<TClass>().cols()),Tidx>::value
                            &&  std::is_convertible<decltype(std::declval<TClass>().matrix()),Mat<Tval>>::value
                            &&  std::is_invocable<decltype(&TClass::getMatrix),TClass,Mat<Tval>&>::value ;
} ;


// TODO: We might at some point want to do the same for low-rank matrices


// Some static assertions on crossables that should satisfy the crossable concept
//  might remove later on

using value_t = double ;
using complex_t = std::complex<value_t> ;

using V1m1 = FnSp::LinDiscFnSpace<value_t> ;
using Slp = Krnl::SlpKernel<value_t> ;
using BilForm = Blf::BilinearForm<V1m1,Slp> ;

static_assert(is_crossable<LR::LRMatrix<complex_t>, complex_t>::value && "LR::LRMatrix<complex_t> is not crossable") ;
static_assert(is_crossable<LR::SVDMatrix<complex_t>,complex_t>::value && "LR::SVDMatrix<complex_t> is not crossable") ;

static_assert(is_crossable<Blf::BlfCrossable<BilForm,value_t,true>, complex_t>::value && "Blf::BlfCrossable<BilForm,true> is not crossable") ;
static_assert(is_crossable<Blf::BlfCrossable<BilForm,complex_t,false>,complex_t>::value && "Blf::BlfCrossable<BilForm,false> is not crossable") ;

static_assert(is_crossable<Cross::SumLRCrossable<value_t>,value_t>::value && "Cross::SumLRCrossable<value_t> is not crossable") ;


} // namespace Util

#endif // _TYPE_TRAITS_HH