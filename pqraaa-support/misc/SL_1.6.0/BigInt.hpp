/*
 * Copyright (c) 2007-2014, A. N. Yzelman,   Utrecht University 2007-2011;
 *                                                    KU Leuven 2011-2014.
 *                          R. H. Bisseling, Utrecht University 2007-2014.
 * 
 * This file is part of the Sparse Library.
 * 
 * This library was developed under supervision of Prof. dr. Rob H. Bisseling at
 * Utrecht University, from 2007 until 2011. From 2011-2014, development continued 
 * at KU Leuven, where Prof. dr. Dirk Roose contributed significantly to the ideas 
 * behind the newer parts of the library code.
 * 
 *     The Sparse Library is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by the
 *     Free Software Foundation, either version 3 of the License, or (at your
 *     option) any later version.
 * 
 *     The Sparse Library is distributed in the hope that it will be useful, but
 *     WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 *     or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 *     for more details.
 * 
 *     You should have received a copy of the GNU General Public License along
 *     with the Sparse Library. If not, see <http://www.gnu.org/licenses/>.
 */


/* @file
 * File created by:
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2012.
 */


#ifndef _H_BIG_INT
#define _H_BIG_INT

/** A 128-bit integer, with overloaded comparison operators. */
struct BigInt {

	/** Most significant bits */
	unsigned long int high;

	/** Least significant bits. */
	unsigned long int low;

	/** Default constructor (sets fields to zero). */
	BigInt(): high( 0 ), low( 0 ) {}

	/**
	 * Direct constructor.
	 * 
	 * @param a Most significant bits.
	 * @param b Least significant bits.
	 */
	BigInt( unsigned long int a, unsigned long int b ): high( a ), low( b ) {}

	/** Assignment operator. @param other the integer this integer should equal. */
	void operator=( const BigInt &other ) {
		high = other.high;
		low  = other.low;
	}

	/** @return unsinged long int cast of this integer. */
	operator unsigned long int() {
		if( high == 0ul )
			return low;
		else
			return high;
	}

	/** @return unsinged int cast of this integer. */
	operator unsigned int() {
		if( high == 0ul )
			return static_cast< unsigned int >( low );
		else
			return static_cast< unsigned int >( high );
	}

	/** @return unsigned short int cast of this integer. */
	operator unsigned short int() {
		if( high == 0ul )
			return static_cast< unsigned short int >( low );
		else
			return static_cast< unsigned short int >( high );
	}

	/** @return unsigned char cast of this integer. */
	operator unsigned char() {
		if( high == 0ul )
			return static_cast< unsigned char >( low );
		else
			return static_cast< unsigned char >( high );
	}

};

/**
 * Compare function used for quicksort on an array of bigints
 * 
 * @param a The pointer to the integer which will be compared.
 * @param b The pointer to the integer to compare with.
 * 
 * @return 0 if equal, 1 if a > b, -1 if b > a.
 */
int bigint_compare( const void *a, const void *b ) {
	const BigInt left  = *(BigInt*)a;
	const BigInt right = *(BigInt*)b;
	if( left.high == right.high ) {
		if( left.low == right.low ) {
			return 0;
		} else {
			return (left.low > right.low) ? 1 : -1;
		}
	} else {
		return (left.high > right.high) ? 1 : -1;
	}
}

/**
 * Comparison operator.
 *
 * @param left Left integer in the comparison.
 * @param right Right integer in the comparison.
 * @return whether left > right.
 */
bool operator>( const BigInt &left, const BigInt &right ) {
	return bigint_compare( &left, &right ) > 0;
}

/**
 * Comparison operator.
 *
 * @param left Left integer in the comparison.
 * @param right Right integer in the comparison.
 * @return wheter left < right.
 */
bool operator<( const BigInt &left, const BigInt &right ) {
	return bigint_compare( &left, &right ) < 0;
}

/**
 * Comparison operator.
 *
 * @param left Left integer in the comparison.
 * @param right Right integer in the comparison.
 * @return whether left == right.
 */
bool operator==( const BigInt &left, const BigInt &right ) {
	return left.high == right.high && left.low == right.low;
}

/**
 * Comparison operator.
 *
 * @param left Left integer in the comparison.
 * @param right Right integer in the comparison.
 * @return whether left != right.
 */
bool operator!=( const BigInt &left, const BigInt &right ) {
	return left.high != right.high || left.low != right.low;
}

/**
 * Comparison operator.
 *
 * @param left Left integer in the comparison.
 * @param right Right integer in the comparison.
 * @return whether left >= right.
 */
bool operator>=( const BigInt &left, const BigInt &right ) {
	return left == right || ( left > right );
}

/**
 * Comparison operator.
 *
 * @param left Left integer in the comparison.
 * @param right Right integer in the comparison.
 * @return whether left <= right.
 */
bool operator<=( const BigInt &left, const BigInt &right ) {
	return left == right || ( left < right );
}

#endif

