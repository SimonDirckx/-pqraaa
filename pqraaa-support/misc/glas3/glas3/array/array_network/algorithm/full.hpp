//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_full_hpp
#define glas3_array_array_network_algorithm_full_hpp

#include <glas3/concept/is.hpp>

#include <glas3/array/array_network/concept/tensor_network.hpp>

#include <glas3/array/type/vector_wrapper.hpp>
#include <glas3/array/dense_array/type/shape_index.hpp>
#include <glas3/array/dense_array/type/range.hpp>
#include <glas3/array/dense_array/container/dense_array.hpp>
#include <glas3/array/dense_array/type/block_selection_index.hpp>

#include <type_traits>
#include <vector>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <tuple>

namespace glas3 {

template <typename X>
typename std::enable_if< is<TensorNetwork, X>::value, dense_array<typename X::value_type> >::type
full( X const& x ) {

	auto nodes2arrays = x.nodes2arrays() ;
	std::map<typename X::ndims_type, boost::shared_ptr<dense_array<typename X::value_type>>> new_nodes2arrays ;

	for ( auto contraction_data: x.contraction_sequence() ) {

		// construct new uninitialized node_array
		boost::shared_ptr<dense_array<typename X::value_type>> new_node_array = boost::make_shared<dense_array<typename X::value_type>>( no_init(), contraction_data.new_node_shape ) ;

		// perform contraction: extensive implementation because of loop pruning
		if ( contraction_data.nodes2contract.size() == 1 ) {
			auto node = *contraction_data.nodes2contract.begin() ;
			if ( new_nodes2arrays.count( node ) ) {
				for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
					typename X::value_type sum = 0 ;
					for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in() ) {
						sum += (*new_nodes2arrays.at( node ))( *contraction_data.nodes2indices[node] ) ;
					}
					(*new_node_array)(*contraction_data.outer_index) = sum ;
				}
			}
			else {
				for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
					typename X::value_type sum = 0 ;
					for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in() ) {
						sum += (*nodes2arrays.at( node ))( *contraction_data.nodes2indices[node] ) ;
					}
					(*new_node_array)(*contraction_data.outer_index) = sum ;
				}
			}
		}
		else {
			auto node0 = *contraction_data.nodes2contract.begin() ;
			auto node1 = *(++contraction_data.nodes2contract.begin()) ;
			dense_vector<typename X::value_type> product( no_init(), contraction_data.index->lin_size() ) ;
			if ( new_nodes2arrays.count( node0 ) ) {
				if ( new_nodes2arrays.count( node1 ) ) {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*new_nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*new_nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
				else {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*new_nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
			}
			else {
				if ( new_nodes2arrays.count( node1 ) ) {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*new_nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
				else {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
			}
		}

		// update nodes2arrays
		for ( auto node: contraction_data.nodes2contract ) {
			if ( new_nodes2arrays.count( node) ) {
				new_nodes2arrays.erase( node ) ;
			}
			else {
				nodes2arrays.erase( node ) ;
			}
		}
		new_nodes2arrays[contraction_data.new_node] = new_node_array ;
	}

	return new_nodes2arrays.begin()->second->shallow_copy() ;
}

template <typename X, typename J>
typename std::enable_if< is<TensorNetwork, X>::value && is<Array, J>::value && std::is_same< typename X::array_type, dense_array<typename X::value_type> >::value, typename X::value_type >::type
full_entry( X const& x, J const& ind ) {

	auto nodes2arrays = x.nodes2arrays() ;
	std::map<typename X::ndims_type, boost::shared_ptr<dense_array<typename X::value_type>>> new_nodes2arrays ;
	(*x.contraction_sequence_entry().front().entry_index) = ind ;

	for ( auto contraction_data: x.contraction_sequence_entry() ) {

		// construct new uninitialized node_array
		boost::shared_ptr<dense_array<typename X::value_type>> new_node_array = boost::make_shared<dense_array<typename X::value_type>>( no_init(), contraction_data.new_node_shape ) ;

		// perform contraction: extensive implementation because of loop pruning
		if ( contraction_data.nodes2contract.size() == 1 ) {
			auto node = *contraction_data.nodes2contract.begin() ;
			if ( new_nodes2arrays.count( node ) ) {
				for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
					typename X::value_type sum = 0 ;
					for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in() ) {
						sum += (*new_nodes2arrays.at( node ))( *contraction_data.nodes2indices[node] ) ;
					}
					(*new_node_array)(*contraction_data.outer_index) = sum ;
				}
			}
			else {
				for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
					typename X::value_type sum = 0 ;
					for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in() ) {
						sum += (*nodes2arrays.at( node ))( *contraction_data.nodes2indices[node] ) ;
					}
					(*new_node_array)(*contraction_data.outer_index) = sum ;
				}
			}
		}
		else {
			auto node0 = *contraction_data.nodes2contract.begin() ;
			auto node1 = *(++contraction_data.nodes2contract.begin()) ;
			dense_vector<typename X::value_type> product( no_init(), contraction_data.index->lin_size() ) ;
			if ( new_nodes2arrays.count( node0 ) ) {
				if ( new_nodes2arrays.count( node1 ) ) {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*new_nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*new_nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
				else {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*new_nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
			}
			else {
				if ( new_nodes2arrays.count( node1 ) ) {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*new_nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
				else {
					for ( contraction_data.outer_index->reset(); contraction_data.outer_index->overflow_count() < 1; contraction_data.outer_index->inc_lin_in() ) {
						typename X::ndims_type i = 0 ;
						for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
							product[i] = (*nodes2arrays.at( node0 ))( *contraction_data.nodes2indices[node0] )
									   * (*nodes2arrays.at( node1 ))( *contraction_data.nodes2indices[node1] )  ;
						}
						for( auto node_it = ++(++contraction_data.nodes2contract.begin()); node_it != contraction_data.nodes2contract.end(); ++node_it ) {
							if ( new_nodes2arrays.count( *node_it ) ) {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*new_nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
							else {
								i = 0 ;
								for ( contraction_data.index->reset(); contraction_data.index->overflow_count() < 1; contraction_data.index->inc_lin_in(), ++i ) {
									product[i] *= (*nodes2arrays.at( *node_it ))( *contraction_data.nodes2indices[*node_it] ) ;
								}
							}
						}
						typename X::value_type sum = 0 ;
						for ( i = 0; i < product.size(); ++i ) {
							sum += product[i] ;
						}
						(*new_node_array)(*contraction_data.outer_index) = sum ;
					}
				}
			}
		}

		// update nodes2arrays
		for ( auto node: contraction_data.nodes2contract ) {
			if ( new_nodes2arrays.count( node) ) {
				new_nodes2arrays.erase( node ) ;
			}
			else {
				nodes2arrays.erase( node ) ;
			}
		}
		new_nodes2arrays[contraction_data.new_node] = new_node_array ;
	}

	return (*new_nodes2arrays.begin()->second)[0] ;
}

} // namespace glas3

#endif

