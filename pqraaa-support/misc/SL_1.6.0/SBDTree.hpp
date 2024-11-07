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


/*
 * File created by:
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2010.
 */


#ifndef _H_SBDTREE
#define _H_SBDTREE

#include <cstdlib>
#include <vector>
#include <iostream>
#include <limits.h>
#include <assert.h>

/** Models a Separated Block Diagonal tree structure. */
class SBDTree {

	protected:
		
		/** Array s.t. parent[i] returns the parent ID of block i. */
		unsigned long int *parent;

		/** Array s.t. left_child[i] returns the left child ID of block i. */
		unsigned long int *left_child;

		/** Array s.t. right_child[i] returns the right child ID of block i. */
		unsigned long int *right_child;

		/** Array s.t. r_lo[i] returns the lower row bound of block i. */
		unsigned long int *r_lo;

		/** Array s.t. r_hi[i] returns the upper row bound of block i. */
		unsigned long int *r_hi;

		/** Array s.t. c_lo[i] returns the lower column bound of block i. */
		unsigned long int *c_lo;

		/** Array s.t. c_hi[i] returns the upper column bound of block i. */
		unsigned long int *c_hi;

		/** Which node ID corresponds to the root. */
		unsigned long int root;

		/** The total number of tree nodes. */
		unsigned long int nodes;

		/** Whether the root node is set. */
		char root_is_set;

		/** Integer corresponding to non-existing nodes. */
		static const unsigned long int NO_SUCH_NODE = ULONG_MAX;

		/**
		 * Builds the SBD tree using three input vectors.
		 * 
		 * @param hierarchy The SBD hierarchy vector.
		 * @param r_bounds  The row bounds of each SBD block.
		 * @param c_bounds  The column bounds of each SBD block.
		 */
		void build( std::vector< unsigned long int > &hierarchy,
				std::vector< unsigned long int > &r_bounds,
				std::vector< unsigned long int > &c_bounds );
	public:

		/** Base constructor */
		SBDTree( std::vector< unsigned long int > &r_hierarchy, std::vector< unsigned long int > &c_hierarchy,
				std::vector< unsigned long int > &r_bounds,
				std::vector< unsigned long int > &c_bounds );

		/** Base constructor. Warning: avoids some assertions! */
		SBDTree( std::vector< unsigned long int > &hierarchy,
				std::vector< unsigned long int > &r_bounds,
				std::vector< unsigned long int > &c_bounds );

		/** Base deconstructor. */
		~SBDTree();

		/**
		 * Gets, from a separator node, the bounding box of the nonzeroes contained in the separator.
		 * Note that this is *not* the r_lo/hi,c_lo/hi of the node itself; those are the bounds of the
		 * row-wise and column-wise separators.
		 *
		 * This is a logarithmic operation.
		 *
		 * @param index The separator node ID.
		 * @param r_lo Where to store the lower row bound.
		 * @param r_hi Where to store the upper row bound.
		 * @param c_lo Where to store the lower column bound.
		 * @param c_hi Where to store the upper column bound.
		 */
		void getSeparatorBB( const unsigned long int index,
					unsigned long int &r_lo, unsigned long int &r_hi,
					unsigned long int &c_lo, unsigned long int &c_hi );

		/** Returns the parent of a given node. */
		unsigned long int up( const unsigned long int index );

		/** Returns the left child of a given node. */
		unsigned long int left( const unsigned long int index );

		/** Returns the right child of a given node. */
		unsigned long int right( const unsigned long int index );

		/** Returns the row bounds corresponding to a given node. */
		void rowBounds( const unsigned long int index,
				unsigned long int &r_lo,
				unsigned long int &r_hi );

		/** Returns the column bounds corresponding to a given node. */
		void columnBounds( const unsigned long int index,
					unsigned long int &c_lo,
					unsigned long int &c_hi );

		/** Whether the given node is a leaf node */
		char isLeaf( const unsigned long int index );

		/** Gets the number of nodes */
		unsigned long int size();

		/** Gets the root index */
		unsigned long int getRoot();
};

#endif

