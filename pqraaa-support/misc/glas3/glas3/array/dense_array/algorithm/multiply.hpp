//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_algorithm_multiply_hpp
#define glas3_array_dense_array_algorithm_multiply_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas2/concept/is.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/matrix/concept/dense_matrix.hpp>
#include <glas2/sparse/concept/sparse_matrix.hpp>

#include <glas3/array/dense_array/container/dense_scalar.hpp>
#include <glas3/array/dense_array/container/dense_vector.hpp>
#include <glas3/array/dense_array/container/dense_matrix.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/ttt_1D_operation.hpp>
#include <glas3/array/dense_array/type/range.hpp>
#include <glas3/array/dense_array/type/shape_index.hpp>
#include <glas3/array/dense_array/type/block_selection_index.hpp>
#include <glas3/array/type/vector_wrapper.hpp>
#include <glas3/array/dense_array/type/blocks_array.hpp>

#include <glas3/array/dense_array/algorithm/iostream.hpp>
#include <glas3/array/dense_array/algorithm/linear_index_selection.hpp>

#include <type_traits>
#include <iostream>

namespace glas3 {

template < typename X, typename Y >
typename std::enable_if< is< DenseVector, X >::value && is< DenseVector, Y >::value, dense_scalar<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.size() == y.size() ) ;
	dense_scalar<value_type> r( 0 ) ;
	for ( typename X::size_type i = 0; i < x.size(); ++i ) {
		r[0] += x[i] * y[i] ;
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseMatrix, X >::value && is< DenseVector, Y >::value, dense_vector<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.shape()[1] == y.size() ) ;
	dense_vector<value_type> r( 0, x.shape()[0] ) ;
	typename X::size_type k = 0 ;
	for ( typename Y::size_type j = 0; j < y.size(); ++j ) {
		for ( typename X::size_type i = 0; i < x.shape()[0]; ++i ) {
			r[i] += x[k] * y[j] ;
			++k ;
		}
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseVector, X >::value && is< DenseMatrix, Y >::value, dense_vector<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.size() == y.shape()[0]  ) ;
	dense_vector<value_type> r( 0, y.shape()[1] ) ;
	typename Y::size_type k = 0 ;
	for ( typename Y::size_type j = 0; j < y.shape()[1]; ++j ) {
		for ( typename X::size_type i = 0; i < x.size(); ++i ) {
			r[j] += x[i] * y[k] ;
			++k ;
		}
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseMatrix, X >::value && is< DenseMatrix, Y >::value, dense_matrix<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.shape()[1] == y.shape()[0]  ) ;
	dense_matrix<value_type> r( 0, {x.shape()[0], y.shape()[1]} ) ;
	std::ptrdiff_t ij, ij_0 ;
	typename X::size_type ik ;
	typename Y::size_type kj ;

	for ( typename Y::size_type j = 0; j < y.shape()[1]; ++j ) {
		ij_0 = x.shape()[0] * j ;
		kj = x.shape()[1] * j ;
		for ( typename X::size_type k = 0; k < x.shape()[1]; ++k ) {
			ij = ij_0 ;
			ik = x.shape()[0] * k ;
			for ( typename X::size_type i = 0; i < x.shape()[0]; ++i ) {
				r[ij] += x[ik] * y[kj] ;
				++ij ;
				++ik ;
			}
			++kj ;
		}
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value
                  && ! ( is< DenseVector, X >::value && is< DenseVector, Y >::value )
                  && ! ( is< DenseMatrix, X >::value && is< DenseVector, Y >::value )
                  && ! ( is< DenseVector, X >::value && is< DenseMatrix, Y >::value )
                  && ! ( is< DenseMatrix, X >::value && is< DenseMatrix, Y >::value ), dense_array<typename ttt_1D_operation< X, Y >::value_type> >::type
multiply ( X const& x, Y const& y ) {
	typedef typename ttt_1D_operation< X, Y >::value_type value_type ;

	return dense_array<value_type>( ttt_1D_operation< X, Y >( x, y, x.shape().size() - 1, 0 ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseVector, X >::value && is< DenseVector, Y >::value, dense_scalar<typename ttt_1D_operation< X, Y >::value_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typedef typename ttt_1D_operation< X, Y >::value_type value_type ;

	return dense_scalar<value_type>( ttt_1D_operation< X, Y >( x, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseMatrix, X >::value && is< DenseVector, Y >::value, dense_vector<typename ttt_1D_operation< X, Y >::value_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typedef typename ttt_1D_operation< X, Y >::value_type value_type ;

	return dense_vector<value_type>( ttt_1D_operation< X, Y >( x, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseVector, X >::value && is< DenseMatrix, Y >::value, dense_vector<typename ttt_1D_operation< X, Y >::value_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typedef typename ttt_1D_operation< X, Y >::value_type value_type ;

	return dense_vector<value_type>( ttt_1D_operation< X, Y >( x, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseMatrix, X >::value && is< DenseMatrix, Y >::value, dense_matrix<typename ttt_1D_operation< X, Y >::value_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typedef typename ttt_1D_operation< X, Y >::value_type value_type ;

	return dense_matrix<value_type>( ttt_1D_operation< X, Y >( x, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseArray, X >::value && is< DenseArray, Y >::value
                  && ! ( is< DenseVector, X >::value && is< DenseVector, Y >::value )
                  && ! ( is< DenseMatrix, X >::value && is< DenseVector, Y >::value )
                  && ! ( is< DenseVector, X >::value && is< DenseMatrix, Y >::value )
                  && ! ( is< DenseMatrix, X >::value && is< DenseMatrix, Y >::value ), dense_array<typename ttt_1D_operation< X, Y >::value_type> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typedef typename ttt_1D_operation< X, Y >::value_type value_type ;

	return dense_array<value_type>( ttt_1D_operation< X, Y >( x, y, inner_dim_x, inner_dim_y ) ) ;
}

template < typename X, typename Y >
typename std::enable_if< glas2::is< glas2::SparseMatrix, X >::value && is< DenseVector, Y >::value, dense_vector<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.num_columns() == y.shape()[0]  ) ;
	dense_vector<value_type> r( 0, x.num_rows() ) ;
	for (typename X::size_type i = 0; i < x.num_nz(); ++i) {
		r[x.row(i)] += x.data()(i) * y[x.column(i)] ;
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseVector, X >::value && glas2::is< glas2::SparseMatrix, Y >::value, dense_vector<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.size() == y.num_rows()  ) ;
	dense_vector<value_type> r( 0, y.num_columns() ) ;
	for (typename Y::size_type i = 0; i < y.num_nz(); ++i) {
		r[y.column(i)] += y.data()(i) * x[y.row(i)] ;
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< glas2::is< glas2::SparseMatrix, X >::value && is< DenseVector, Y >::value, dense_vector<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y, typename Y::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	assert( inner_dim_y == 0 ) ;
	assert( inner_dim_x == 0 || inner_dim_x == 1 ) ;

	if ( inner_dim_x == 1 ) { return multiply( x, y ) ; }
	else { return multiply( y, x ) ; }
}

template < typename X, typename Y >
typename std::enable_if< is< DenseVector, X >::value && glas2::is< glas2::SparseMatrix, Y >::value, dense_vector<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename X::ndims_type inner_dim_y ) {
	assert( inner_dim_x == 0 ) ;
	assert( inner_dim_y == 0 || inner_dim_y == 1 ) ;

	if ( inner_dim_y == 0 ) { return multiply( x, y ) ; }
	else { return multiply( y, x ) ; }
}

template < typename X, typename Y >
typename std::enable_if< glas2::is< glas2::SparseMatrix, X >::value && is< DenseMatrix, Y >::value, dense_matrix<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.num_columns() == y.shape()[0]  ) ;
	dense_matrix<value_type> r( 0, {x.num_rows(), y.shape()[1]} ) ;
	for (typename X::size_type i = 0; i < x.num_nz(); ++i) {
		for ( typename Y::size_type j = 0; j < y.shape()[1]; ++j ) {
			r({x.row(i), j}) += x.data()(i) * y({x.column(i), j}) ;
		}
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseMatrix, X >::value && glas2::is< glas2::SparseMatrix, Y >::value, dense_matrix<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.shape()[1] == y.num_rows()  ) ;
	dense_matrix<value_type> r( 0, {x.shape()[0], y.num_columns()} ) ;
	for (typename Y::size_type i = 0; i < y.num_nz(); ++i) {
		for ( typename X::size_type j = 0; j < x.shape()[0]; ++j ) {
			r({j, y.column(i)}) += y.data()(i) * x({j, y.row(i)}) ;
		}
	}
	return r ;
}

template < typename X, typename Y >
typename std::enable_if< glas2::is< glas2::SparseMatrix, X >::value && is< DenseMatrix, Y >::value, dense_matrix<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y, typename Y::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( inner_dim_x == 0 || inner_dim_x == 1 ) ;
	assert( inner_dim_y == 0 || inner_dim_y == 1 ) ;

	if ( inner_dim_x == 1 ) {
		if ( inner_dim_y == 0 ) { return multiply( x, y ) ; }
		else {
			assert( x.num_columns() == y.shape()[1]  ) ;
			dense_matrix<value_type> r( 0, {x.num_rows(), y.shape()[0]} ) ;
			for (typename X::size_type i = 0; i < x.num_nz(); ++i) {
				for ( typename Y::size_type j = 0; j < y.shape()[0]; ++j ) {
					r({x.row(i), j}) += x.data()(i) * y({j, x.column(i)}) ;
				}
			}
			return r ;
		}
	}
	else {
		if ( inner_dim_y == 1 ) { return multiply( y, x ) ; }
		else {
			assert( y.shape()[0] == x.num_rows()  ) ;
			dense_matrix<value_type> r( 0, {y.shape()[1], x.num_columns()} ) ;
			for (typename X::size_type i = 0; i < x.num_nz(); ++i) {
				for ( typename Y::size_type j = 0; j < y.shape()[1]; ++j ) {
					r({j, x.column(i)}) += x.data()(i) * y({x.row(i), j}) ;
				}
			}
			return r ;
		}
	}
}

template < typename X, typename Y >
typename std::enable_if< is< DenseMatrix, X >::value && glas2::is< glas2::SparseMatrix, Y >::value, dense_matrix<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename X::ndims_type inner_dim_y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( inner_dim_x == 0 || inner_dim_x == 1 ) ;
	assert( inner_dim_y == 0 || inner_dim_y == 1 ) ;

	if ( inner_dim_x == 1 ) {
		if ( inner_dim_y == 0 ) { return multiply( x, y ) ; }
		else {
			assert( y.num_columns() == x.shape()[1]  ) ;
			dense_matrix<value_type> r( 0, {x.shape()[0], y.num_rows()} ) ;
			for (typename Y::size_type i = 0; i < y.num_nz(); ++i) {
				for ( typename X::size_type j = 0; j < x.shape()[0]; ++j ) {
					r({j, y.row(i)}) += y.data()(i) * x({j, y.column(i)}) ;
				}
			}
			return r ;
		}
	}
	else {
		if ( inner_dim_y == 1 ) { return multiply( y, x ) ; }
		else {
			assert( x.shape()[0] == y.num_rows()  ) ;
			dense_matrix<value_type> r( 0, {y.num_columns(), x.shape()[1]} ) ;
			for (typename Y::size_type i = 0; i < y.num_nz(); ++i) {
				for ( typename X::size_type j = 0; j < x.shape()[1]; ++j ) {
					r({y.column(i), j}) += y.data()(i) * x({y.row(i), j}) ;
				}
			}
			return r ;
		}
	}
}

template < typename X, typename Y >
typename std::enable_if< glas2::is< glas2::SparseMatrix, X >::value && is< DenseArray, Y >::value && ! is< DenseMatrix, Y >::value && ! is< DenseVector, Y >::value,
                         dense_array<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.num_columns() == y.shape()[0]  ) ;

	std::vector<boost::shared_ptr<vector_wrapper<typename Y::size_type>>> selection_arrays ( y.shape().size() + 1 );
	dense_scalar<typename Y::size_type> i_x( 0 ), i_inner( 0 ) ;

	selection_arrays[0] = boost::make_shared<vector_wrapper<typename Y::size_type>>( i_x ) ;
	selection_arrays[1] = boost::make_shared<vector_wrapper<typename Y::size_type>>( i_inner ) ;
	for ( typename Y::ndims_type k = 1; k < y.shape().size(); ++k ) {
		selection_arrays[k + 1] = boost::make_shared<vector_wrapper<typename Y::size_type>>( range< typename Y::size_type >( 0, y.shape()[k] ) ) ;
	}

	blocks_array<typename Y::size_type> shape_r( {vector_wrapper<typename Y::size_type>( x.num_rows() ), y.shape()[range< typename Y::size_type >( 1, y.shape().size() )]}, {0} ) ;
	blocks_array<typename Y::size_type> shape_index( {vector_wrapper<typename Y::size_type>( x.num_rows() ), vector_wrapper<typename Y::size_type>( y.shape() ) }, {0} ) ;

	block_selection_index<vector_wrapper<typename Y::size_type>> index( shape_index, selection_arrays ) ;
	dense_array<value_type> r( 0, shape_r ) ;

	dense_vector<typename Y::size_type> r_sel( blocks_array<typename Y::size_type> ( {vector_wrapper<typename Y::size_type>( 0 ), range< typename Y::size_type >( 2, y.shape().size() + 1 )}, {0} ) ) ;
	dense_vector<typename Y::size_type> y_sel( range< typename Y::size_type >( 1, y.shape().size() + 1 ) ) ;

	auto j_r = index[r_sel] ;
	auto j_y = index[y_sel] ;

	for ( typename X::size_type i = 0; i < x.num_nz(); ++i ) {
		i_x[0] = x.row(i) ;
		i_inner[0] = x.column(i) ;
		for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
			r(j_r) += x.data()(i) * y(j_y) ;
		}
	}

	return r ;
}

template < typename X, typename Y >
typename std::enable_if< is< DenseArray, X >::value && ! is< DenseMatrix, X >::value && ! is< DenseVector, X >::value && glas2::is< glas2::SparseMatrix, Y >::value,
                         dense_array<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( x.shape()[x.shape().size() - 1] == y.num_rows()  ) ;

	std::vector<boost::shared_ptr<vector_wrapper<typename X::size_type>>> selection_arrays ( x.shape().size() + 1 );
	dense_scalar<typename X::size_type> i_y( 0 ), i_inner( 0 ) ;

	for ( typename X::ndims_type k = 0; k < x.shape().size() - 1; ++k ) {
		selection_arrays[k] = boost::make_shared<vector_wrapper<typename X::size_type>>( range< typename X::size_type >( 0, x.shape()[k] ) ) ;
	}
	selection_arrays[x.shape().size() - 1] = boost::make_shared<vector_wrapper<typename X::size_type>>( i_inner ) ;
	selection_arrays[x.shape().size()] = boost::make_shared<vector_wrapper<typename X::size_type>>( i_y ) ;

	blocks_array<typename X::size_type> shape_r( {x.shape()[range< typename X::size_type >( 0, x.shape().size() - 1 )], vector_wrapper<typename Y::size_type>( y.num_columns() )}, {0} ) ;
	blocks_array<typename X::size_type> shape_index( {vector_wrapper<typename Y::size_type>( x.shape() ), vector_wrapper<typename Y::size_type>( y.num_columns() )}, {0} ) ;

	block_selection_index<vector_wrapper<typename X::size_type>> index( shape_index, selection_arrays ) ;
	dense_array<value_type> r( 0, shape_r ) ;

	dense_vector<typename X::size_type> r_sel( blocks_array<typename X::size_type>( {range< typename X::size_type >( 0, x.shape().size() - 1 ), vector_wrapper<typename X::size_type>( x.shape().size() )}, {0} ) ) ;
	dense_vector<typename X::size_type> x_sel( range< typename X::size_type >( 0, x.shape().size() ) ) ;

	auto j_r = index[r_sel] ;
	auto j_x = index[x_sel] ;

	for ( typename Y::size_type i = 0; i < y.num_nz(); ++i ) {
		i_y[0] = y.column(i) ;
		i_inner[0] = y.row(i) ;
		for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
			r(j_r) += y.data()(i) * x(j_x) ;
		}
	}

	return r ;
}

template < typename X, typename Y >
typename std::enable_if< glas2::is< glas2::SparseMatrix, X >::value && is< DenseArray, Y >::value && ! is< DenseMatrix, Y >::value && ! is< DenseVector, Y >::value,
                         dense_array<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y, typename Y::ndims_type inner_dim_x, typename Y::ndims_type inner_dim_y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( inner_dim_x == 0 || inner_dim_x == 1 ) ;
	assert( inner_dim_y >= 0 && inner_dim_y < y.shape().size() ) ;

	if ( inner_dim_x == 1 ) {
		if ( inner_dim_y == 0 ) { return multiply( x, y ) ; }
		else {
			assert( x.num_columns() == y.shape()[inner_dim_y]  ) ;

			std::vector<boost::shared_ptr<vector_wrapper<typename Y::size_type>>> selection_arrays ( y.shape().size() + 1 );
			dense_scalar<typename Y::size_type> i_x( 0 ), i_inner( 0 ) ;

			selection_arrays[0] = boost::make_shared<vector_wrapper<typename Y::size_type>>( i_x ) ;
			for ( typename Y::ndims_type k = 0; k < y.shape().size(); ++k ) {
				if ( k == inner_dim_y ) {
					selection_arrays[k + 1] = boost::make_shared<vector_wrapper<typename Y::size_type>>( i_inner ) ;
				}
				else {
					selection_arrays[k + 1] = boost::make_shared<vector_wrapper<typename Y::size_type>>( range< typename Y::size_type >( 0, y.shape()[k] ) ) ;
				}
			}

			blocks_array<typename Y::size_type> shape_r( {vector_wrapper<typename Y::size_type>( x.num_rows() ), y.shape()[range< typename Y::size_type >( 0, inner_dim_y )],
				y.shape()[range< typename Y::size_type >( inner_dim_y + 1, y.shape().size() )]}, {0} ) ;
			blocks_array<typename Y::size_type> shape_index( {vector_wrapper<typename Y::size_type>( x.num_rows() ), vector_wrapper<typename Y::size_type>( y.shape() )}, {0} ) ;

			block_selection_index<vector_wrapper<typename Y::size_type>> index( shape_index, selection_arrays ) ;
			dense_array<value_type> r( 0, shape_r ) ;

			dense_vector<typename Y::size_type> r_sel( blocks_array<typename Y::size_type>( {vector_wrapper<typename Y::size_type>( 0 ), range< typename Y::size_type >( 1, inner_dim_y + 1 ),
				range< typename Y::size_type >( inner_dim_y + 2, y.shape().size() + 1 )}, {0} ) ) ;
			dense_vector<typename Y::size_type> y_sel( range< typename Y::size_type >( 1, y.shape().size() + 1 ) ) ;

			auto j_r = index[r_sel] ;
			auto j_y = index[y_sel] ;

			for ( typename X::size_type i = 0; i < x.num_nz(); ++i ) {
				i_x[0] = x.row(i) ;
				i_inner[0] = x.column(i) ;
				for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
					r(j_r) += x.data()(i) * y(j_y) ;
				}
			}

			return r ;
		}
	}
	else {
		if ( inner_dim_y == y.shape().size() - 1 ) { return multiply( y, x ) ; }
		else {
			assert( x.num_rows() == y.shape()[inner_dim_y]  ) ;

			std::vector<boost::shared_ptr<vector_wrapper<typename Y::size_type>>> selection_arrays ( y.shape().size() + 1 );
			dense_scalar<typename Y::size_type> i_x( 0 ), i_inner( 0 ) ;

			selection_arrays[0] = boost::make_shared<vector_wrapper<typename Y::size_type>>( i_x ) ;
			for ( typename Y::ndims_type k = 0; k < y.shape().size(); ++k ) {
				if ( k == inner_dim_y ) {
					selection_arrays[k + 1] = boost::make_shared<vector_wrapper<typename Y::size_type>>( i_inner ) ;
				}
				else {
					selection_arrays[k + 1] = boost::make_shared<vector_wrapper<typename Y::size_type>>( range< typename Y::size_type >( 0, y.shape()[k] ) ) ;
				}
			}

			blocks_array<typename Y::size_type> shape_r( {y.shape()[range< typename Y::size_type >( 0, inner_dim_y )],
				y.shape()[range< typename Y::size_type >( inner_dim_y + 1, y.shape().size() )], vector_wrapper<typename Y::size_type>( x.num_columns() )}, {0} ) ;
			blocks_array<typename Y::size_type> shape_index( {vector_wrapper<typename Y::size_type>( x.num_columns() ), vector_wrapper<typename Y::size_type>( y.shape() )}, {0} ) ;

			block_selection_index<vector_wrapper<typename Y::size_type>> index( shape_index, selection_arrays ) ;
			dense_array<value_type> r( 0, shape_r ) ;

			dense_vector<typename Y::size_type> r_sel( blocks_array<typename Y::size_type>( {range< typename Y::size_type >( 1, inner_dim_y + 1 ),
				range< typename Y::size_type >( inner_dim_y + 2, y.shape().size() + 1 ), vector_wrapper<typename Y::size_type>( 0 )}, {0} ) ) ;
			dense_vector<typename Y::size_type> y_sel( range< typename Y::size_type >( 1, y.shape().size() + 1 ) ) ;

			auto j_r = index[r_sel] ;
			auto j_y = index[y_sel] ;

			for ( typename X::size_type i = 0; i < x.num_nz(); ++i ) {
				i_x[0] = x.column(i) ;
				i_inner[0] = x.row(i) ;
				for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
					r(j_r) += x.data()(i) * y(j_y) ;
				}
			}

			return r ;
		}
	}
}

template < typename X, typename Y >
typename std::enable_if< is< DenseArray, X >::value && ! is< DenseMatrix, X >::value && ! is< DenseVector, X >::value && glas2::is< glas2::SparseMatrix, Y >::value,
                         dense_array<decltype( typename X::value_type() * typename Y::value_type() )> >::type
multiply ( X const& x, Y const& y, typename X::ndims_type inner_dim_x, typename X::ndims_type inner_dim_y ) {
	typedef decltype( typename X::value_type() * typename Y::value_type() ) value_type ;

	assert( inner_dim_x >= 0 && inner_dim_x < x.shape().size() ) ;
	assert( inner_dim_y == 0 || inner_dim_y == 1 ) ;

	if ( inner_dim_y == 1 ) {
		if ( inner_dim_x == 0 ) { return multiply( y, x ) ; }
		else {
			assert( y.num_columns() == x.shape()[inner_dim_x]  ) ;

			std::vector<boost::shared_ptr<vector_wrapper<typename X::size_type>>> selection_arrays ( x.shape().size() + 1 );
			dense_scalar<typename X::size_type> i_y( 0 ), i_inner( 0 ) ;

			for ( typename X::ndims_type k = 0; k < x.shape().size(); ++k ) {
				if ( k == inner_dim_x ) {
					selection_arrays[k] = boost::make_shared<vector_wrapper<typename X::size_type>>( i_inner ) ;
				}
				else {
					selection_arrays[k] = boost::make_shared<vector_wrapper<typename X::size_type>>( range< typename X::size_type >( 0, x.shape()[k] ) ) ;
				}
			}
			selection_arrays[x.shape().size()] = boost::make_shared<vector_wrapper<typename X::size_type>>( i_y ) ;

			blocks_array<typename X::size_type> shape_r( {x.shape()[range< typename X::size_type >( 0, inner_dim_x )], vector_wrapper<typename X::size_type>( y.num_rows() ),
				x.shape()[range< typename X::size_type >( inner_dim_x + 1, x.shape().size() )]}, {0} ) ;
			blocks_array<typename X::size_type> shape_index( {vector_wrapper<typename X::size_type>( x.shape() ), vector_wrapper<typename X::size_type>( y.num_rows() )}, {0} ) ;

			block_selection_index<vector_wrapper<typename X::size_type>> index( shape_index, selection_arrays ) ;
			dense_array<value_type> r( 0, shape_r ) ;

			dense_vector<typename X::size_type> r_sel( blocks_array<typename X::size_type>( {range< typename X::size_type >( 0, inner_dim_x ),
				vector_wrapper<typename X::size_type>( x.shape().size() ), range< typename X::size_type >( inner_dim_x + 1, x.shape().size() )}, {0} ) ) ;
			dense_vector<typename X::size_type> x_sel( range< typename X::size_type >( 0, x.shape().size() ) ) ;

			auto j_r = index[r_sel] ;
			auto j_x = index[x_sel] ;

			for ( typename Y::size_type i = 0; i < y.num_nz(); ++i ) {
				i_y[0] = y.row(i) ;
				i_inner[0] = y.column(i) ;
				for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
					r(j_r) += y.data()(i) * x(j_x) ;
				}
			}

			return r ;
		}
	}
	else {
		if ( inner_dim_x == x.shape().size() - 1 ) { return multiply( x, y ) ; }
		else {
			assert( y.num_rows() == x.shape()[inner_dim_x]  ) ;

			std::vector<boost::shared_ptr<vector_wrapper<typename X::size_type>>> selection_arrays ( x.shape().size() + 1 );
			dense_scalar<typename X::size_type> i_y( 0 ), i_inner( 0 ) ;

			for ( typename X::ndims_type k = 0; k < x.shape().size(); ++k ) {
				if ( k == inner_dim_x ) {
					selection_arrays[k] = boost::make_shared<vector_wrapper<typename X::size_type>>( i_inner ) ;
				}
				else {
					selection_arrays[k] = boost::make_shared<vector_wrapper<typename X::size_type>>( range< typename X::size_type >( 0, x.shape()[k] ) ) ;
				}
			}
			selection_arrays[x.shape().size()] = boost::make_shared<vector_wrapper<typename X::size_type>>( i_y ) ;

			blocks_array<typename X::size_type> shape_r( {x.shape()[range< typename X::size_type >( 0, inner_dim_x )], vector_wrapper<typename X::size_type>( y.num_columns() ),
				x.shape()[range< typename X::size_type >( inner_dim_x + 1, x.shape().size() )]}, {0} ) ;
			blocks_array<typename X::size_type> shape_index( {vector_wrapper<typename X::size_type>( x.shape() ), vector_wrapper<typename X::size_type>( y.num_columns() )}, {0} ) ;

			block_selection_index<vector_wrapper<typename X::size_type>> index( shape_index, selection_arrays ) ;
			dense_array<value_type> r( 0, shape_r ) ;

			dense_vector<typename X::size_type> r_sel( blocks_array<typename X::size_type>( {range< typename X::size_type >( 0, inner_dim_x ),
				vector_wrapper<typename X::size_type>( x.shape().size() ), range< typename X::size_type >( inner_dim_x + 1, x.shape().size() )}, {0} ) ) ;
			dense_vector<typename X::size_type> x_sel( range< typename X::size_type >( 0, x.shape().size() ) ) ;

			auto j_r = index[r_sel] ;
			auto j_x = index[x_sel] ;

			for ( typename Y::size_type i = 0; i < y.num_nz(); ++i ) {
				i_y[0] = y.column(i) ;
				i_inner[0] = y.row(i) ;
				for ( index.reset(); index.overflow_count() < 1; index.inc_lin_in() ) {
					r(j_r) += y.data()(i) * x(j_x) ;
				}
			}

			return r ;
		}
	}
}

} // namespace glas3

#endif
