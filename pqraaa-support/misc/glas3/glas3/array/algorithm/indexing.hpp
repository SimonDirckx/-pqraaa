//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_algorithm_indexing_hpp
#define glas3_array_algorithm_indexing_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>
#include <glas3/array/dense_array/concept/contiguous_dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <cassert>
#include <type_traits>
#include <iostream>
#include <map>
#include <set>

namespace glas3 {

template <typename list_type>
typename list_type::size_type bisect_right(list_type const& a, typename list_type::value_type x, typename list_type::size_type lo, typename list_type::size_type hi) {
    /*Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e <= x, and all e in
    a[i:] have e > x.  So if x already appears in the list, a.insert(x) will
    insert just after the rightmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    */

	typename list_type::size_type mid ;

    while ( lo < hi ) {
        mid = (lo + hi) / 2 ;
        if (x < a[mid]) { hi = mid ; }
        else { lo = mid + 1 ; }
    }
    return lo ;
}

template <typename list_type>
typename list_type::size_type bisect_right( list_type const& a, typename list_type::value_type x ) {
	return bisect_right( a, x, 0, a.size() ) ;
}

template <typename list_type>
typename list_type::size_type bisect_left(list_type const& a, typename list_type::value_type x, typename list_type::size_type lo, typename list_type::size_type hi) {

	/* Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
	*/

	typename list_type::size_type mid ;

    while ( lo < hi ) {
        mid = (lo + hi) / 2 ;
        if (a[mid] < x) { lo = mid + 1 ; }
        else { hi = mid ; }
    }
    return lo ;
}

template <typename list_type>
typename list_type::size_type bisect_left( list_type const& a, typename list_type::value_type x ) {
	return bisect_left( a, x, 0, a.size() ) ;
}

//template <typename shape_type>
//typename shape_type::value_type
//multi_index2index( shape_type const& shape ) {
//    return 0 ;
//}
//
//template <typename shape_type, typename J, typename... I>
//typename std::enable_if< std::is_integral<J>::value, typename shape_type::value_type >::type
//multi_index2index( shape_type const& shape, J const& j, I const&... i ) {
//	typename shape_type::value_type size_current_dim =  shape[shape.size() - sizeof...(I) - 1] ;
//    assert( j >= 0 && j < size_current_dim ) ;
//
//	return j + size_current_dim * multi_index2index( shape, i...) ;
//}

//template <typename shape_type, typename dim_factor_type>
//typename dim_factor_type::value_type
//multi_index2index_dim_factor( shape_type const& shape, dim_factor_type const& dim_factor) {
//    return 0 ;
//}
//
//template <typename shape_type, typename dim_factor_type, typename J, typename... I>
//typename std::enable_if< std::is_integral<J>::value, typename dim_factor_type::value_type >::type
//multi_index2index_dim_factor( shape_type const& shape, dim_factor_type const& dim_factor, J const& j, I const&... i ) {
//	typename shape_type::value_type      size_current_dim   =  shape[shape.size() - sizeof...(I) - 1] ;
//	typename dim_factor_type::value_type factor_current_dim =  dim_factor[shape.size() - sizeof...(I) - 1] ;
//    assert( j >= 0 && j < size_current_dim ) ;
//
//	return factor_current_dim * j + multi_index2index_dim_factor( shape, dim_factor, i...) ;
//}

template <typename shape_type, typename index_type >
typename std::enable_if< is<Array, shape_type>::value && is<Array, index_type>::value, typename shape_type::value_type >::type
multi_index2index_initializer_list( shape_type const& shape, index_type const& j ) {

	typename shape_type::value_type i = 0 ;
	auto i_j = j.size() ;
	if ( i_j != 0 ) {
		--i_j ;
		auto k = shape.size() - 1 ;

		assert( j[i_j] >= 0 && j[i_j] < shape[k] ) ;
		i = j[i_j] ;

		--i_j ; --k ;
		for ( ; k >= 0; --i_j, --k ) {
			assert( j[i_j] >= 0 && j[i_j] < shape[k] ) ;
			i = i * shape[k] + j[i_j] ;
		}
	}

	return i ;
}

template <typename shape_target_type, typename shape_type, typename index_type >
typename std::enable_if< is<Array, shape_target_type>::value && is<Array, shape_type>::value && is<Array, index_type>::value, typename shape_type::value_type >::type
multi_index2index_mod_shape_initializer_list( shape_target_type const& shape_target, shape_type const& shape, index_type const& j ) {

	typename shape_type::value_type i = 0 ;
	auto i_j = j.size() ;
	if ( i_j != 0 ) {
		--i_j ;
		auto k = shape.size() - 1 ;

		assert( j[i_j] >= 0 && j[i_j] < shape[k] ) ;
		i = j[i_j] % shape_target[k] ;

		--i_j ; --k ;
		for ( ; k >= 0; --i_j, --k ) {
			assert( j[i_j] >= 0 && j[i_j] < shape[k] ) ;
			i = i * shape_target[k] + j[i_j] % shape_target[k] ;
		}
	}

	return i ;
}

template <typename shape_type, typename selection_arrays_type, typename index_type >
typename std::enable_if< is<Array, shape_type>::value && is<Array, index_type>::value, typename shape_type::value_type >::type
multi_index2index_selection_arrays_initializer_list( shape_type const& shape, shape_type const& shape_selection, selection_arrays_type const& selection_arrays, index_type const& j ) {

	typename shape_type::value_type i = 0 ;
	auto i_j = j.size() ;
	if ( i_j != 0 ) {
		--i_j ;
		auto k = shape_selection.size() - 1 ;

		assert( j[i_j] >= 0 && j[i_j] < shape_selection[k] ) ;
		typename shape_type::value_type dummy = (*(selection_arrays[k]))[j[i_j]] ;
		assert( dummy >= 0 && dummy < shape[k] ) ;
		i = dummy ;

		--i_j ; --k ;
		for ( ; k >= 0; --i_j, --k ) {
			assert( j[i_j] >= 0 && j[i_j] < shape_selection[k] ) ;
			dummy = (*(selection_arrays[k]))[j[i_j]] ;
			assert( dummy >= 0 && dummy < shape[k] ) ;
			i = i * shape[k] + dummy ;
		}
	}

	return i ;
}

template <typename shape_type, typename dim_factor_type, typename index_type >
typename std::enable_if< is<Array, shape_type>::value && is<Array, dim_factor_type>::value && is<Array, index_type>::value, typename shape_type::value_type >::type
multi_index2index_dim_factor_initializer_list( shape_type const& shape, dim_factor_type const& dim_factor, index_type const& j ) {

	typename shape_type::value_type i = 0 ;

	if ( j.size() != 0 ) {
		typename index_type::ndims_type i_j = 0 ;
		for ( typename shape_type::size_type k = 0; k < shape.size(); ++i_j, ++k ) {
			assert( j[i_j] >= 0 && j[i_j] < shape[k] ) ;
			i += j[i_j] * dim_factor[k] ;
		}
	}

	return i ;
}

template < typename dim_factor_type, typename shape_type, typename shape_selection_type, typename selection_arrays_type, typename index_type >
typename std::enable_if< is<Array, shape_type>::value && is<Array, shape_selection_type>::value && is<Array, dim_factor_type>::value && is<Array, index_type>::value, typename shape_type::value_type >::type
multi_index2index_dim_factor_selection_arrays_initializer_list( shape_type const& shape, dim_factor_type const& dim_factor, shape_selection_type const& shape_selection, selection_arrays_type const& selection_arrays, index_type const& j ) {

    typename shape_type::value_type i = 0, dummy ;

    if ( j.size() != 0 ) {
    	typename index_type::ndims_type i_j = 0 ;
    	for ( typename shape_type::size_type k = 0; k < shape_selection.size(); ++i_j, ++k ) {
    		assert( j[i_j] >= 0 && j[i_j] < shape_selection[k] ) ;
    		dummy = (*(selection_arrays[k]))[j[i_j]] ;
    		assert( dummy >= 0 && dummy < shape[k] ) ;
    		i += dummy * dim_factor[k] ;
    	}
    }

	return i ;
}

template < typename D, typename I, typename M >
typename std::enable_if< is<Array, D>::value && std::is_integral<I>::value >::type
index2multi_index( D const& dim_factor_in, I i_in, M& j_in ) {

    std::ldiv_t dv{};
    dv.rem = i_in ;
    typename D::value_type i_out = 0 ;

    if ( dim_factor_in.size() != 0 ) {
    	typename D::size_type k = bisect_right( dim_factor_in, i_in + 1 ), l ;
    	for ( typename M::size_type ii = j_in.size() - 1; ii >= k; --ii ) { j_in[ii] = 0 ; }
    	while ( k > 0 ) {
    		l = k - 1 ;
    		if ( dv.rem >= dim_factor_in[l] ) {
    			dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );
    			j_in[l] = dv.quot ;
    		}
    		else {
    			j_in[l] = 0 ;
    		}
    		--k ;
    	}
    }
}

template < typename D, typename I >
typename std::enable_if< is<Array, D>::value && std::is_integral<I>::value, typename D::value_type >::type
index2multi_index2index( D const& dim_factor_target, D const& dim_factor, I i ) {

    std::ldiv_t dv{};
    dv.rem = i ;
    typename D::value_type j = 0 ;

    if ( dim_factor.size() != 0 ) {
    	typename D::size_type k = bisect_right( dim_factor, i + 1 ), l ;
    	while ( k > 0 ) {
    		l = k - 1 ;
    		if ( dv.rem >= dim_factor[l] ) {
    			dv = std::ldiv( dv.rem, long(dim_factor[l]) );
    			j += dv.quot * dim_factor_target[l] ;
    		}
    		--k ;
    	}
    }

    return j ;
}

template < typename D, typename I, typename M >
typename std::enable_if< is<Array, D>::value && std::is_integral<I>::value, typename D::value_type >::type
index2multi_index2index( D const& dim_factor_out, D const& dim_factor_in, I i_in, M& j_in ) {

    std::ldiv_t dv{};
    dv.rem = i_in ;
    typename D::value_type i_out = 0 ;

    if ( dim_factor_in.size() != 0 ) {
    	typename D::size_type k = bisect_right( dim_factor_in, i_in + 1 ), l ;
    	for ( typename M::size_type ii = j_in.size() - 1; ii >= k; --ii ) { j_in[ii] = 0 ; }
    	while ( k > 0 ) {
    		l = k - 1 ;
    		if ( dv.rem >= dim_factor_in[l] ) {
    			dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );
    			j_in[l] = dv.quot ;
    			i_out += dv.quot * dim_factor_out[l] ;
    		}
    		else {
    			j_in[l] = 0 ;
    		}
    		--k ;
    	}
    }

    return i_out ;
}

template < typename D, typename I, typename M, typename shape_out_type >
typename std::enable_if< is<Array, D>::value && std::is_integral<I>::value, typename D::value_type >::type
index2multi_index2index_mod_shape( D const& dim_factor_out, shape_out_type const& shape_out, D const& dim_factor_in, I i_in, M& j_in, M& j_out ) {

    std::ldiv_t dv{};
    dv.rem = i_in ;
    typename D::value_type i_out = 0 ;

    if ( dim_factor_in.size() != 0 ) {
    	typename D::size_type k = bisect_right( dim_factor_in, i_in + 1 ), l ;
    	for ( typename M::size_type ii = j_in.size() - 1; ii >= k; --ii ) { j_in[ii] = 0 ; j_out[ii] = 0 ; }
    	while ( k > 0 ) {
    		l = k - 1 ;
    		if ( dv.rem >= dim_factor_in[l] ) {
    			dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );
    			j_in[l] = dv.quot ;
    			j_out[l] = j_in[l] % shape_out[l] ;
    			i_out += j_out[l] * dim_factor_out[l] ;
    		}
    		else {
    			j_in[l] = 0 ;
    			j_out[l] = 0 ;
    		}
    		--k ;
    	}
    }

    return i_out ;
}

template < typename S, typename D, typename A, typename I, typename M >
typename std::enable_if< is<Array, S>::value && is<Array, D>::value && std::is_integral<I>::value, typename S::value_type >::type
index2multi_index2index_selection_arrays( S const& shape_out, D const& dim_factor_in, A const& selection_arrays, I i_in, M& j_in, M& j_out ) {

    std::ldiv_t dv{};
    dv.rem = i_in ;
    typename D::size_type l ;
    typename S::value_type i_out = 0 ;

    for ( typename D::size_type k = dim_factor_in.size(); k > 0; --k ) {
       	l = k - 1 ;
       	if ( dv.rem >= dim_factor_in[l] ) {
       		dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );
       	}
       	j_in[l] = dv.quot ;
       	auto dummy = (*(selection_arrays[l]))[dv.quot] ;
       	j_out[l] = dummy ;
       	assert( dummy >= 0 && dummy < shape_out[l] ) ;
       	dv.quot = 0 ;

       	i_out = dummy + i_out * shape_out[l] ;
    }

    return i_out ;
}

template < typename S, typename D, typename A, typename I >
typename std::enable_if< is<Array, S>::value && is<Array, D>::value && std::is_integral<I>::value, typename S::value_type >::type
index2multi_index2index_selection_arrays( S const& shape_out, D const& dim_factor_in, A const& selection_arrays, I i_in ) {

    std::ldiv_t dv{};
    dv.rem = i_in ;
    typename D::size_type l ;
    typename S::value_type i_out = 0 ;

    for ( typename D::size_type k = dim_factor_in.size(); k > 0; --k ) {
       	l = k - 1 ;
       	if ( dv.rem >= dim_factor_in[l] ) {
       		dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );
       	}
       	auto dummy = (*(selection_arrays[l]))[dv.quot] ;
       	assert( dummy >= 0 && dummy < shape_out[l] ) ;
       	dv.quot = 0 ;

       	if ( k < dim_factor_in.size() ) { i_out = dummy + i_out * shape_out[l] ; }
       	else { i_out = dummy ; }
    }

    return i_out ;
}

template <typename D, typename S, typename M, typename I1, typename I2 >
void index2multi_index_memory( D const& dim_factor, S const& shape, M& j, I1 const& i_in_new, I2& i_in ) {
	std::map<typename S::size_type, typename M::value_type> changed_entries ;
	typename D::size_type k, l ;

	if ( dim_factor.size() != 0 ) {
		if ( i_in_new != i_in ) {
			if ( i_in_new > i_in ) {
				std::ldiv_t dv{} ;
				dv.rem = i_in_new - i_in ;
				k = bisect_right( dim_factor, dv.rem + 1 );
				while ( k > 0 ) {
					l = k - 1;
					if ( dv.rem >= dim_factor[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor[l]) );

						while ( dv.quot > 0 ) {
							j[l] += dv.quot ;
							if ( j[l] >= shape[l] ) {

								j[l] -= shape[l] ;

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			else {
				std::ldiv_t dv{} ;
				dv.rem = i_in - i_in_new ;
				k = bisect_right( dim_factor, dv.rem + 1 ) ;
				while ( k > 0 ) {
					l = k - 1 ;
					if ( dv.rem >= dim_factor[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor[l]) );

						while ( dv.quot > 0 ) {
							if ( j[l] < dv.quot ) {

								j[l] += shape[l] - dv.quot;

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								j[l] -= dv.quot ;
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			i_in = i_in_new ;
		}
	}
}

template <typename D, typename S, typename M, typename I1, typename I2 >
void index2multi_index2index_memory( D const& dim_factor_target, D const& dim_factor, S const& shape, M& j, I1 i_in_new, I2& i_in, I2& i_out ) {
	std::map<typename S::size_type, typename M::value_type> changed_entries ;
	typename D::size_type k, l ;

	if ( dim_factor.size() != 0 ) {
		if ( i_in_new != i_in ) {
			if ( i_in_new > i_in ) {
				std::ldiv_t dv{} ;
				dv.rem = i_in_new - i_in ;
				k = bisect_right( dim_factor, dv.rem + 1 );
				while ( k > 0 ) {
					l = k - 1;
					if ( dv.rem >= dim_factor[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j[l] ) ) ;

						while ( dv.quot > 0 ) {
							j[l] += dv.quot ;
							if ( j[l] >= shape[l] ) {

								j[l] -= shape[l] ;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			else {
				std::ldiv_t dv{} ;
				dv.rem = i_in - i_in_new ;
				k = bisect_right( dim_factor, dv.rem + 1 ) ;
				while ( k > 0 ) {
					l = k - 1 ;
					if ( dv.rem >= dim_factor[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j[l] ) ) ;

						while ( dv.quot > 0 ) {
							if ( j[l] < dv.quot ) {

								j[l] += shape[l] - dv.quot;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								j[l] -= dv.quot ;
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			for (auto e = changed_entries.begin(); e != changed_entries.end(); ++e) {
				auto dummy = j[e->first] ;
				if ( dummy > e->second ) { i_out += ( dummy - e->second ) * dim_factor_target[e->first] ; }
				else if ( dummy < e->second ) { i_out -= ( e->second - dummy ) * dim_factor_target[e->first] ; }
			}
			i_in = i_in_new ;
		}
	}
}

template <typename D, typename S, typename A, typename M, typename I1, typename I2 >
void index2multi_index_selection_arrays_memory( S const& shape_target, D const& dim_factor, S const& shape, A const& selection_arrays, M& j_in, M& j_out, I1 i_in_new, I2& i_in) {
	std::set<typename S::size_type> changed_entries ;
	typename D::size_type k, l ;

	if ( dim_factor.size() != 0 ) {
		if ( i_in_new != i_in ) {
			if ( i_in_new > i_in ) {
				std::ldiv_t dv{} ;
				dv.rem = i_in_new - i_in ;
				k = bisect_right( dim_factor, dv.rem + 1 );
				while ( k > 0 ) {
					l = k - 1;
					if ( dv.rem >= dim_factor[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor[l]) ) ;

						changed_entries.insert( l ) ;

						while ( dv.quot > 0 ) {
							j_in[l] += dv.quot ;
							if ( j_in[l] >= shape[l] ) {

								j_in[l] -= shape[l] ;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( l + 1 ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			else {
				std::ldiv_t dv{} ;
				dv.rem = i_in - i_in_new ;
				k = bisect_right( dim_factor, dv.rem + 1 ) ;
				while ( k > 0 ) {
					l = k - 1 ;
					if ( dv.rem >= dim_factor[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor[l]) );

						changed_entries.insert( l ) ;

						while ( dv.quot > 0 ) {
							if ( j_in[l] < dv.quot ) {

								j_in[l] += shape[l] - dv.quot;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( l + 1 ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								j_in[l] -= dv.quot ;
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			for (auto e = changed_entries.begin(); e != changed_entries.end(); ++e) {
				j_out[*e] = (*(selection_arrays[*e]))[j_in[*e]] ;
				assert( j_out[*e] >= 0 && j_out[*e] < shape_target[*e] ) ;
			}
			i_in = i_in_new ;
		}
	}
}

template <typename D, typename S, typename A, typename M, typename I1, typename I2 >
void index2multi_index2index_selection_arrays_memory(D const& dim_factor_out, S const& shape_out, D const& dim_factor_in, S const& shape_in, A const& selection_arrays, M& j_in, I1 i_in_new, I2& i_in, I2& i_out) {
	std::map<typename S::size_type, typename M::value_type> changed_entries ;
	typename D::size_type k, l ;

	if ( dim_factor_in.size() != 0 ) {
		if ( i_in_new != i_in ) {
			if ( i_in_new > i_in ) {
				std::ldiv_t dv{} ;
				dv.rem = i_in_new - i_in ;
				k = bisect_right( dim_factor_in, dv.rem + 1 );
				while ( k > 0 ) {
					l = k - 1;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							j_in[l] += dv.quot ;
							if ( j_in[l] >= shape_in[l] ) {

								j_in[l] -= shape_in[l] ;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			else {
				std::ldiv_t dv{} ;
				dv.rem = i_in - i_in_new ;
				k = bisect_right( dim_factor_in, dv.rem + 1 ) ;
				while ( k > 0 ) {
					l = k - 1 ;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							if ( j_in[l] < dv.quot ) {

								j_in[l] += shape_in[l] - dv.quot;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								j_in[l] -= dv.quot ;
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			for (auto e = changed_entries.begin(); e != changed_entries.end(); ++e){
				auto dummy1 = (*(selection_arrays[e->first]))[j_in[e->first]] ;
				auto dummy2 = (*(selection_arrays[e->first]))[e->second] ;
				assert( dummy1 >= 0 && dummy1 < shape_out[e->first] ) ;
				if ( dummy1 > dummy2 ) { i_out += ( dummy1 - dummy2 ) * dim_factor_out[e->first] ; }
				else if ( dummy1 < dummy2 ) { i_out -= ( dummy2 - dummy1 ) * dim_factor_out[e->first] ; }
			}
			i_in = i_in_new ;
		}
	}
}

template <typename D, typename S, typename A, typename M, typename I1, typename I2 >
void index2multi_index2index_selection_arrays_memory(D const& dim_factor_out, S const& shape_out, D const& dim_factor_in, S const& shape_in, A const& selection_arrays, M& j_in, M& j_out, I1 i_in_new, I2& i_in, I2& i_out) {
	std::map<typename S::size_type, typename M::value_type> changed_entries ;
	typename D::size_type k, l ;

	if ( dim_factor_in.size() != 0 ) {
		if ( i_in_new != i_in ) {
			if ( i_in_new > i_in ) {
				std::ldiv_t dv{} ;
				dv.rem = i_in_new - i_in ;
				k = bisect_right( dim_factor_in, dv.rem + 1 );
				while ( k > 0 ) {
					l = k - 1;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							j_in[l] += dv.quot ;
							if ( j_in[l] >= shape_in[l] ) {

								j_in[l] -= shape_in[l] ;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			else {
				std::ldiv_t dv{} ;
				dv.rem = i_in - i_in_new ;
				k = bisect_right( dim_factor_in, dv.rem + 1 ) ;
				while ( k > 0 ) {
					l = k - 1 ;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							if ( j_in[l] < dv.quot ) {

								j_in[l] += shape_in[l] - dv.quot;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								j_in[l] -= dv.quot ;
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			for (auto e = changed_entries.begin(); e != changed_entries.end(); ++e){
				auto dummy1 = (*(selection_arrays[e->first]))[j_in[e->first]] ;
				j_out[e->first] = dummy1 ;
				auto dummy2 = (*(selection_arrays[e->first]))[e->second] ;
				assert( dummy1 >= 0 && dummy1 < shape_out[e->first] ) ;
				if ( dummy1 > dummy2 ) { i_out += ( dummy1 - dummy2 ) * dim_factor_out[e->first] ; }
				else if ( dummy1 < dummy2 ) { i_out -= ( dummy2 - dummy1 ) * dim_factor_out[e->first] ; }
			}
			i_in = i_in_new ;
		}
	}
}

template <typename D, typename S, typename shape_target_type, typename M, typename I1, typename I2 >
void index2multi_index2index_mod_shape_memory(D const& dim_factor_out, shape_target_type const& shape_out, D const& dim_factor_in, S const& shape_in, M& j_in, I1 i_in_new, I2& i_in, I2& i_out) {
	std::map<typename S::size_type, typename M::value_type> changed_entries ;
	typename D::size_type k, l ;

	if ( dim_factor_in.size() != 0 ) {
		if ( i_in_new != i_in ) {
			if ( i_in_new > i_in ) {
				std::ldiv_t dv{} ;
				dv.rem = i_in_new - i_in ;
				k = bisect_right( dim_factor_in, dv.rem + 1 );
				while ( k > 0 ) {
					l = k - 1;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							j_in[l] += dv.quot ;
							if ( j_in[l] >= shape_in[l] ) {

								j_in[l] -= shape_in[l] ;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			else {
				std::ldiv_t dv{} ;
				dv.rem = i_in - i_in_new ;
				k = bisect_right( dim_factor_in, dv.rem + 1 ) ;
				while ( k > 0 ) {
					l = k - 1 ;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							if ( j_in[l] < dv.quot ) {

								j_in[l] += shape_in[l] - dv.quot;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								j_in[l] -= dv.quot ;
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			for (auto e = changed_entries.begin(); e != changed_entries.end(); ++e) {
				auto dummy = j_in[e->first] % shape_out[e->first] ;
				auto dummy2 = e->second % shape_out[e->first] ;
				if ( dummy > dummy2 ) { i_out += ( dummy - dummy2 ) * dim_factor_out[e->first] ; }
				else if ( dummy < dummy2 ) { i_out -= ( dummy2 - dummy ) * dim_factor_out[e->first] ; }
			}
			i_in = i_in_new ;
		}
	}
}

template <typename D, typename S, typename shape_target_type, typename M, typename I1, typename I2 >
void index2multi_index2index_mod_shape_memory(D const& dim_factor_out, shape_target_type const& shape_out, D const& dim_factor_in, S const& shape_in, M& j_in, M& j_out, I1 i_in_new, I2& i_in, I2& i_out) {
	std::map<typename S::size_type, typename M::value_type> changed_entries ;
	typename D::size_type k, l ;

	if ( dim_factor_in.size() != 0 ) {
		if ( i_in_new != i_in ) {
			if ( i_in_new > i_in ) {
				std::ldiv_t dv{} ;
				dv.rem = i_in_new - i_in ;
				k = bisect_right( dim_factor_in, dv.rem + 1 );
				while ( k > 0 ) {
					l = k - 1;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							j_in[l] += dv.quot ;
							if ( j_in[l] >= shape_in[l] ) {

								j_in[l] -= shape_in[l] ;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			else {
				std::ldiv_t dv{} ;
				dv.rem = i_in - i_in_new ;
				k = bisect_right( dim_factor_in, dv.rem + 1 ) ;
				while ( k > 0 ) {
					l = k - 1 ;
					if ( dv.rem >= dim_factor_in[l] ) {
						dv = std::ldiv( dv.rem, long(dim_factor_in[l]) );

						changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l, j_in[l] ) ) ;

						while ( dv.quot > 0 ) {
							if ( j_in[l] < dv.quot ) {

								j_in[l] += shape_in[l] - dv.quot;

								if ( changed_entries.count(l + 1) == 0 ) {
									changed_entries.insert( std::pair<typename S::size_type, typename M::value_type>( l + 1, j_in[l + 1] ) ) ;
								}

								dv.quot = 1 ;
								l += 1 ;
							}
							else {
								j_in[l] -= dv.quot ;
								dv.quot = 0 ;
							}
						}
					}
					--k ;
				}
			}
			for (auto e = changed_entries.begin(); e != changed_entries.end(); ++e) {
				auto dummy = j_in[e->first] % shape_out[e->first] ;
				j_out[e->first] = dummy ;
				auto dummy2 = e->second % shape_out[e->first] ;
				if ( dummy > dummy2 ) { i_out += ( dummy - dummy2 ) * dim_factor_out[e->first] ; }
				else if ( dummy < dummy2 ) { i_out -= ( dummy2 - dummy ) * dim_factor_out[e->first] ; }
			}
			i_in = i_in_new ;
		}
	}
}

} // namespace glas3

#endif
