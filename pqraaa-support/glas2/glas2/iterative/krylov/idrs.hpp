//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas_toolbox_iterative_krylov_idrs_hpp
#define glas_toolbox_iterative_krylov_idrs_hpp

#include <glas/config.hpp>
#undef GLAS_INCLUDE_LEVEL
#define GLAS_MY_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX
#include <glas/concept/check_include_level.hpp>

#include <glas/glas.hpp>

#include <glas/operators/all.hpp>
#include <glas/container/dense_vector.hpp>
#include <glas/container/dense_matrix.hpp>
#include <glas/toolbox/iterative/krylov/tolerance.hpp>
#include <glas/iterator/begin.hpp>
#include <glas/iterator/end.hpp>
#include <glas/algorithm/abs.hpp>
#include <glas/algorithm/norm_2.hpp>
#include <glas/algorithm/norm_fro.hpp>
#include <glas/algorithm/column.hpp>
#include <glas/algorithm/vector_range.hpp>
#include <glas/type/random_seed.hpp>
#include <glas/algorithm/row_range.hpp>
//#include <glas/algorithm/matrix_range.hpp>
#include <glas/algorithm/column_range.hpp>
#include <glas/algorithm/row.hpp>
#include <glas/algorithm/conj.hpp>
#include <glas/algorithm/herm.hpp>
#include <glas/algorithm/upper.hpp>
#include <glas/algorithm/diagonal.hpp>
#include <glas/concept/value_type.hpp>
#include <glas/concept/integral_equal.hpp>
#include <cassert>
#include <algorithm>
#include <cmath>

#include <iostream>

#undef GLAS_INCLUDE_LEVEL
#define GLAS_INCLUDE_LEVEL GLAS_INCLUDE_TOOLBOX

namespace glas {

struct idrs_options: options {
	unsigned int n_vectors_;

	inline idrs_options() :
		n_vectors_(4) {
	}
};

// Reference implementation of
// An elegant IDR(s) variant that efficiently exploits bi-orthogonality
// properties, M. van Gijzen and P. Sonneveld, Tech. Rep., Delft University of
// Technology.

// Op and Prec are here separate operators
template<typename X, typename Y, typename Prec, typename Op, typename Report>
void idrs(Op const& op, Prec const& prec, X& x, Y const& b, Report& report,
		idrs_options const& opt) {
	BOOST_STATIC_ASSERT(
		(boost::is_same<
			matrix_vector_product_tag,
			typename linear_operator_category<Op>::type
		>::value)
	) ;

        using glas::operators::operator- ;
        using glas::operators::operator* ;
        using glas::operators::operator+ ;

	typedef typename value_type<X>::type value_type;
        typedef typename abs_result_type<value_type>::type real_value_type ;
	typedef dense_vector<value_type> vector_type;
	typedef dense_matrix<value_type> basis_type;

	const int s = opt.n_vectors_;

	// Allocate memory for matrices.
	basis_type Q(size(x), s);
	basis_type G(size(x), s);
	basis_type U(size(x), s);
	basis_type M( eye<value_type,int,int>(s,s) );
	G = U = 0.0;

	// Allocate memory for vectors.
	vector_type res(size(x));
	vector_type v(size(x));
	vector_type t( size(x) );
	vector_type c(s);
	vector_type f(s);
	vector_type alpha(s);

	// Select a random s-dimensional basis and orthogonalize it using stabilized
	// Gram-Schmidt.
	random_seed<real_value_type> seed;
        randomize( Q, seed ) ;
	/*for (int i = 0; i < Q.num_rows() * Q.num_columns(); ++i) {
		seed.var_gen(*(Q.storage_ptr() + i));
	}*/
	for (int j = 0; j < s; ++j) {
		c[range(0,j)] = herm( Q[all()][j] ) * Q[all()][range(0,j)];
		Q[all()][j] -= Q[all()][range(0,j)] * c[range(0,j)];
		Q[all()][j] /= norm_2(Q[all()][j]);
	}

	// Compute initial residual.
	op(x, res);
	res = b - res;
	const real_value_type res_norm_0 = norm_2(res);
	report( res_norm_0, 0 );
	if( res_norm_0 == 0.0 ) {
		return;
	}

	// Solve the system iteratively.
	value_type omega = 1;
	for (unsigned int it = 0; it <= opt.max_mat_vec_ / (s+1); ++it) {

		// Compute s independent column vectors for G_it
		f = herm(Q) * res;
		for (int k = 0; k < s; ++k) {

			// Compute v \in G_it \cap S
			c[range(k,s)] =
				inverse( lower(M[range(k,s)][range(k,s)]) ) * f[range(k,s)];
			v  = res - (G[all()][range(k,s)] * c[range(k,s)]);
			prec( v );

			U[all()][k]  = U[all()][k] * c(k);
			U[all()][k] += (U[all()][range(k+1,s)] * c[range(k+1,s)]);
			U[all()][k] += (omega*v);
			/*
			 * NOTE Following statement is not a correct equivalent for above
			 * statements due to lazy evaluation:
			 * 		U[all][k] = (U[all()][range(k,s)] * c[range(k,s)]) +
			 *   		(omega * v)
			 */

			// Compute g_k \in G_it
			op(U[all()][k], t);
			G[all()][k] = t;

			// Bi-orthogonalize g_k and u_k w.r.t. Q
			value_type alpha;
			for(int j = 0; j < k; ++j) {
				alpha =	( herm(Q[all()][j])*G[all()][k] ) / M(j,j);
				G[all()][k] -= alpha * G[all()][j];
				U[all()][k] -= alpha * U[all()][j];
			}

			// Update M
			M[ range(k,s) ][k] = herm( Q[all()][range(k,s)] ) * G[all()][k];

			// Bi-orthogonalize residual w.r.t. Q
			const value_type beta = f(k) / M(k,k);
			res -= beta * G[all()][k];
			x   += beta * U[all()][k];

			const real_value_type res_norm = norm_2(res);
			report(res_norm, it*(s+1) + k + 1);

			if( res_norm <= tolerance( opt, res_norm_0 ) ) {
				return;
			}

			// Update f.
			if(k < s-1) {
				f[ range(k+1,s) ] -= beta * M[range(k+1,s)][k];
			}
		}

		// Update residual.
		v = res;
		prec( v );
		op(v, t);

		const real_value_type ntmp = norm_2(t);
		const value_type dptr = dot(herm(t), res);
		const real_value_type rho = glas::abs( dptr / ( ntmp*norm_2(res) ) );
		omega = dptr / (ntmp*ntmp);
		if( rho < 0.7 ) {
			omega *= 0.7/rho;
		}

		res -= omega * t;
		x   += omega * v;

		real_value_type res_norm = norm_2(res);
		report(res_norm, (it+1)*(s+1) );

		// Stopping criterion.
                if( res_norm <= tolerance( opt, res_norm_0 ) ) {
			return;
		}
	}
} // idrs

} // namespace glas

#endif
