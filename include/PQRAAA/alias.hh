#ifndef _ALIAS_HH
#define _ALIAS_HH
//
// Project     : Alias
// File        : alias.hh
// Description : Collects all aliases used throughout the package
// Author      : Kobe Bruyninckx
// Copyright   : KU Leuven Dept. CS 2023-2024. 
//

#include <vector>
#include <fstream>
#include <cstring>
#include <exception>
#include <iostream>
// #include <filesystem> // New since C++17
#include <sys/types.h> // Necessary for macOS/OS X

#include <BEACHpack/Core>
#include <Eigen/SparseCore>

// TO-DO: Maybe we can do multiple alias files?
namespace Alias {

// Eigen dynamic matrices
template< typename Tval >
using Mat = Eigen::Matrix<Tval,Eigen::Dynamic,Eigen::Dynamic>;
template< typename Tval >
using Mat2X = Eigen::Matrix<Tval,2,Eigen::Dynamic>;
template< typename Tval >
using MatX2 = Eigen::Matrix<Tval,Eigen::Dynamic,2>;
template< typename Tval >
using Mat3X = Eigen::Matrix<Tval,3,Eigen::Dynamic>;
template< typename Tval >
using MatX3 = Eigen::Matrix<Tval,Eigen::Dynamic,3>;
template< typename Tval >
using Mat4X = Eigen::Matrix<Tval,4,Eigen::Dynamic>;
template< typename Tval >
using MatX4 = Eigen::Matrix<Tval,Eigen::Dynamic,4>;

// Eigen fixed-size matrices
template< typename Tval >
using Mat22 = Eigen::Matrix<Tval,2,2>;
template< typename Tval >
using Mat24 = Eigen::Matrix<Tval,2,4>;
template< typename Tval >
using Mat32 = Eigen::Matrix<Tval,3,2>;
template< typename Tval >
using Mat33 = Eigen::Matrix<Tval,3,3>;
template< typename Tval >
using Mat34 = Eigen::Matrix<Tval,3,4>;
template< typename Tval >
using Mat44 = Eigen::Matrix<Tval,4,4>;

// Eigen vectors
template< typename Tval >
using Vec = Eigen::Matrix<Tval,Eigen::Dynamic,1>;
template< typename Tval >
using Vec2 = Eigen::Matrix<Tval,2,1>;
template< typename Tval >
using Vec3 = Eigen::Matrix<Tval,3,1>;
template< typename Tval >
using Vec4 = Eigen::Matrix<Tval,4,1>;

// Eigen sparse matrices and related types
template< typename Tval >
using spMat = Eigen::SparseMatrix<Tval>;
template< typename Tval >
using Triplet = Eigen::Triplet<Tval>;
using c_int = uint ; // integer type for compact sparse matrices

// Eigen indexing types etc.
using Tidx = Eigen::Index ;  // This is just long instead of unsigned long (size_t)
} // namespace Alias
#endif // _ALIAS_HH
