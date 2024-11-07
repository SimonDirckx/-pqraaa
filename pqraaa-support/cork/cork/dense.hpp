//  (C) Copyright Karl Meerbergen & Dries De Samblanx, 2021.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_dense_hpp
#define cork_dense_hpp

#ifdef cork_nonlinear_matrix_hpp
#error #include <cork/dense.hpp> should come before #include <cork/nonlinear_matrix.hpp>
#endif

#ifdef cork_polynomial_matrix_hpp
#error #include <cork/dense.hhp> should be included before #include <cork/polynomial_matrix.hpp>
#endif

#include <cork/coefficient_matrices/dense.hpp>
#include <cork/coefficient_matrices/any_glas.hpp>
#include <cork/linear_solver/any_glas.hpp>
#include <cork/linear_solver/dense.hpp>

#endif
