
#ifndef _H_MCBSP_OPS
#define _H_MCBSP_OPS

#include <limits>
#include <cstdlib>

namespace mcbsp {

/**
 * Namespace containing various operations.
 *
 * It contains standard implementations of the following operations:
 *     +, -, *, /, min, max.
 * All these operations are defined for all standard numeric value types.
 *
 * In a mathematical sense, these operations correspond to additive actions as
 * used for groups, but also includes multiplicative actions as used for (semi-)
 * rings.
 *
 * Each operation defines the `apply(i,j)'-function, which applies the action
 * \f$ ij \f$ and overwrites \f$ i \f$ with the result of that action. For
 * performance reasons, each operation also defines `apply(x,y,n)' where x, y
 * are arrays of length n. The array x will be overwritten with the result of
 * mapping the action on each pair of elements in x and y, i.e.,
 *     \f$ \forall i \in \{0,1,\ldots,n\}, x_i := x_iy_i. \f$
 * Each operation defines whether it is additive and/or multiplicative.
 *
 * For each additive action, the identity element is defined. For each
 * multiplicative action, the zero element is defined. Let \f$ S \f$ be the
 * set of all elements the given action is defined on. Then:
 * 
 *     -a identity element \f$ e \in S \f$ is a unique element for which
 *      \f$ \forall i \in S, ie=i. \f$
 *     -a zero element \f$ z \in S \f$ is the unique element for which
 *      \f$ \forall i \in S, iz=z. \f$
 *
 * 
 */
namespace ops {

	/**
	 * Standard additions of numbers.
	 *
	 * Addition is an additive action, but not a multiplicative one since there
	 * exists no zero element. Addition under numbers forms a group.
	 *
	 * @tparam ValueType which type is used to represent numbers.
	 */
	template< typename ValueType >
	struct Add {

		/** Addition on numbers forms a group. */
		static constexpr bool additive = true;

		/** Addition on numbers as multiplicative action cannot (no identity). */
		static constexpr bool multiplicative = false;

		/** Identity under the addition `+' group action, i.e., 0. */
		static constexpr ValueType identity = static_cast< ValueType >(0);

		/** Given \f$ i,j \f$, replaces \f$ i \f$ by \f$ ij \f$. */
		static void apply( ValueType &left, const ValueType &right ) {
			apply( &left, &right, 1 );
		}

		/** 
		 * Given \f$ n \in \mathbbm{N} \f$ and \f$ x,y \in S^n\f$, replaces each
		 * \f$ x_i \f$ by \f$ x_iy_i \f$, \f$ \forall i \in \{0,1,\ldots,n-1\}. \f$
		 */
		static void apply( ValueType * __restrict__ left, const ValueType * __restrict__ right, const size_t n ) {
			for( size_t i = 0; i < n; ++i ) {
				left[ i ] += right[ i ];
			}
		}
	};

	/**
	 * Standard subtraction of numbers.
	 *
	 * Subtraction is an additive action, but not a multiplicative one since there
	 * exists no zero element. Subtraction under numbers forms a group. It is the
	 * inverse action of standard addition.
	 *
	 * @tparam ValueType which type is used to represent numbers.
	 */
	template< typename ValueType >
	struct Sub {

		/** Subtraction on numbers forms a group. */
		static constexpr bool additive = true;

		/** Subtraction on numbers as multiplicative action cannot (no identity). */
		static constexpr bool multiplicative = false;

		/** Identity under the subtraction `-' group action, i.e., 0. */
		static constexpr ValueType identity = static_cast< ValueType >(0);

		/** Given \f$ i,j \f$, replaces \f$ i \f$ by \f$ ij \f$. */
		static void apply( ValueType &left, const ValueType &right ) {
			apply( &left, &right, 1 );
		}

		/** 
		 * Given \f$ n \in \mathbbm{N} \f$ and \f$ x,y \in S^n\f$, replaces each
		 * \f$ x_i \f$ by \f$ x_iy_i \f$, \f$ \forall i \in \{0,1,\ldots,n-1\}. \f$
		 */
		static void apply( ValueType * __restrict__ left, const ValueType * __restrict__ right, const size_t n ) {
			for( size_t i = 0; i < n; ++i ) {
				left[ i ] -= right[ i ];
			}
		}
	};

	/**
	 * Standard multiplication of numbers.
	 *
	 * Multiplication is an additive action and also a multiplicative action.
	 * Multiplication under numbers forms a ring.
	 *
	 * @tparam ValueType which type is used to represent numbers.
	 */
	template< typename ValueType >
	struct Mul {

		/** Multiplication on numbers forms a group. */
		static constexpr bool additive = true;

		/** Multiplication on numbers as multiplicative action can form a ring. */
		static constexpr bool multiplicative = true;

		/** The zero value under standard multiplication. */
		static constexpr ValueType zero     = static_cast< ValueType >(0);

		/** Identity under the multiplication group action, i.e., 1. */
		static constexpr ValueType identity = static_cast< ValueType >(1);

		/** Given \f$ i,j \f$, replaces \f$ i \f$ by \f$ ij \f$. */
		static void apply( ValueType &left, const ValueType &right ) {
			apply( &left, &right, 1 );
		}

		/** 
		 * Given \f$ n \in \mathbbm{N} \f$ and \f$ x,y \in S^n\f$, replaces each
		 * \f$ x_i \f$ by \f$ x_iy_i \f$, \f$ \forall i \in \{0,1,\ldots,n-1\}. \f$
		 */
		static void apply( ValueType * __restrict__ left, const ValueType * __restrict__ right, const size_t n ) {
			for( size_t i = 0; i < n; ++i ) {
				left[ i ] *= right[ i ];
			}
		}
	};

	/**
	 * Standard division of numbers.
	 *
	 * Division is not associative, since \f$ (6/3)/2 = 1 \neq 6/(3/2)=12/3=4. \f$
	 * Thus division on numbers does not form a group. When excluding 0 from the
	 * set of numbers, division still defines an identity element \f$ e=1 \f$ s.t.
	 * \f$ xe = x \f$ for any number \f$ x \neq 0 \f$. The action then must remain
	 * right-sided, since division also does not commute (\f$ a/b \neq b/a \f$).
	 *
	 * If division is instead applied left-sided, and if again 0 is excluded from
	 * the set of numbers the left-sided division is applied on, then division
	 * defines a zero element \f$ z=0 \f$ s.t. \f$ zx = z \f$ for any number
	 * \f$ x \neq 0. \f$
	 *
	 * Hence even though division under numbers does not yield a group or (semi-)
	 * ring, this operation still defines a zero and identity element.
	 *
	 * @tparam ValueType which type is used to represent numbers.
	 */
	template< typename ValueType >
	struct Div {

		/** Division on numbers does not form a group. */
		static constexpr bool additive = false;

		/** Division on numbers as a multiplicative action cannot form a ring. */
		static constexpr bool multiplicative = false;

		/** The zero value under left-sided division, i.e., 0. */
		static constexpr ValueType zero     = static_cast< ValueType >(0);

		/** Identity under right-sided division, i.e., 1. */
		static constexpr ValueType identity = static_cast< ValueType >(1);

		/** Given \f$ i,j \f$, replaces \f$ i \f$ by \f$ ij \f$. */
		static void apply( ValueType &left, const ValueType &right ) {
			apply( &left, &right, 1 );
		}

		/** 
		 * Given \f$ n \in \mathbbm{N} \f$ and \f$ x,y \in S^n\f$, replaces each
		 * \f$ x_i \f$ by \f$ x_iy_i \f$, \f$ \forall i \in \{0,1,\ldots,n-1\}. \f$
		 */
		static void apply( ValueType * __restrict__ left, const ValueType * __restrict__ right, const size_t n ) {
			for( size_t i = 0; i < n; ++i ) {
				left[ i ] /= right[ i ];
			}
		}
	};

	/**
	 * Taking the maximum of two numbers.
	 *
	 * Taking the maximum is an additive action and also a multiplicative one.
	 * It also forms a ring on numbers. On two numbers \f$ i, j \in S \f$, it
	 * is defined as \f$ ij = \max\{i,j\}. \f$
	 *
	 * @tparam ValueType which type is used to represent numbers.
	 */
	template< typename ValueType >
	struct Max {

		/**
		 * Taking the maximum is an additive action and forms a group under the set
		 * of numbers.
		 */
		static constexpr bool additive = true;

		/**
		 * Taking the maximum is a multiplicative action and forms a group under the
		 * set of numbers.
		 */
		static constexpr bool multiplicative = true;

		/** The zero value under the maximum action. */
		static constexpr ValueType zero     = std::numeric_limits< ValueType >::max();

		/** Identity under the maximum group action. */
		static constexpr ValueType identity = std::numeric_limits< ValueType >::min();

		/** Given \f$ i,j \f$, replaces \f$ i \f$ by \f$ ij \f$. */
		static void apply( ValueType &left, const ValueType &right ) {
			apply( &left, &right, 1 );
		}

		/** 
		 * Given \f$ n \in \mathbbm{N} \f$ and \f$ x,y \in S^n\f$, replaces each
		 * \f$ x_i \f$ by \f$ x_iy_i \f$, \f$ \forall i \in \{0,1,\ldots,n-1\}. \f$
		 */
		static void apply( ValueType * __restrict__ left, const ValueType * __restrict__ right, const size_t n ) {
			for( size_t i = 0; i < n; ++i ) {
				if( left[ i ] < right[ i ] ) {
					left[ i ] = right[ i ];
				}
			}
		}
	};

	/**
	 * Taking the minimum of two numbers.
	 *
	 * Taking the minimum is an additive action and also a multiplicative one.
	 * It also forms a ring on numbers, and is the inverse action of taking
	 * the maximum. For any two numbers \f$ i, j \in S \f$, the minimum action
	 * is defined as \f$ ij = \min \{i,j\}. \f$
	 *
	 * @tparam ValueType which type is used to represent numbers.
	 */
	template< typename ValueType >
	struct Min {

		/**
		 * Taking the minimum is an additive action and forms a group under the set
		 * of numbers.
		 */
		static constexpr bool additive = true;

		/**
		 * Taking the minimum is a multiplicative action and forms a group under the
		 * set of numbers.
		 */
		static constexpr bool multiplicative = true;

		/** The zero value under the minimum action. */
		static constexpr ValueType zero     = std::numeric_limits< ValueType >::min();

		/** Identity under the minimum group action. */
		static constexpr ValueType identity = std::numeric_limits< ValueType >::max();

		/** Given \f$ i,j \f$, replaces \f$ i \f$ by \f$ ij \f$. */
		static void apply( ValueType &left, const ValueType &right ) {
			apply( &left, &right, 1 );
		}

		/** 
		 * Given \f$ n \in \mathbbm{N} \f$ and \f$ x,y \in S^n\f$, replaces each
		 * \f$ x_i \f$ by \f$ x_iy_i \f$, \f$ \forall i \in \{0,1,\ldots,n-1\}. \f$
		 */
		static void apply( ValueType * __restrict__ left, const ValueType * __restrict__ right, const size_t n ) {
			for( size_t i = 0; i < n; ++i ) {
				if( left[ i ] > right[ i ] ) {
					left[ i ] = right[ i ];
				}
			}
		}
	};

} //ops
} //mcbsp

#endif

