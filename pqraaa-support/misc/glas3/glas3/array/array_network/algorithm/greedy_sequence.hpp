//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_algorithm_greedy_sequence_hpp
#define glas3_array_array_network_algorithm_greedy_sequence_hpp

#include <glas3/array/dense_array/type/range.hpp>
#include <glas3/array/array_network/type/nodes_contraction.hpp>

#include <type_traits>
#include <vector>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <tuple>

namespace glas3 {

template <typename NE, typename EN, typename ES, typename IE, typename OE >
typename std::deque<nodes_contraction<typename IE::value_type, typename ES::mapped_type>> greedy_sequence( NE nodes2edges, EN edges2nodes,
		ES edges2sizes, IE inner_edges, OE outer_edges ) {

	typedef typename IE::value_type ndims_type ;
	typedef typename ES::mapped_type size_type ;

	ndims_type i, j, k ;

	std::deque<nodes_contraction<ndims_type, size_type>> contraction_sequence ;

	std::set<ndims_type> set_inner_edges( inner_edges.begin(), inner_edges.end() ) ;
	std::set<ndims_type> set_outer_edges( outer_edges.begin(), outer_edges.end() ) ;

	std::map<ndims_type, std::set<ndims_type>>         current_edges2nodes ;
	std::map<ndims_type, std::vector<ndims_type>>      current_nodes2edges ;
	for ( auto edges2nodes_it: edges2nodes ) {
		auto edge = edges2nodes_it.first ;
		for ( auto nodes_it: edges2nodes_it.second ) {
			current_edges2nodes[edge].insert( std::get<0>( nodes_it ) ) ;
		}
	}
	for ( auto nodes2edges_it: nodes2edges ) {
		auto node = nodes2edges_it.first ;
		for ( auto edge: nodes2edges_it.second ) {
			current_nodes2edges[node].push_back( edge ) ;
		}
	}

	while ( set_inner_edges.size() > 0 || current_nodes2edges.size() > 1 ) {
		std::map<ndims_type, nodes_contraction<ndims_type, size_type>> contractions ;
		if ( set_inner_edges.size() == 0 ) {
			// store in local variable
			contractions[0].new_node = current_nodes2edges.begin()->first ;
			for ( auto node2edges: current_nodes2edges ) {
				contractions.at( 0 ).nodes2contract.insert( node2edges.first ) ;
				if ( node2edges.first < contractions.at( 0 ).new_node ) { contractions.at( 0 ).new_node = node2edges.first ; }
			}
			for ( auto node: contractions.at( 0 ).nodes2contract ) {
				for ( auto edge: current_nodes2edges.at( node ) ) {
					contractions.at( 0 ).edges2keep.insert( edge ) ;
				}
			}
		}
		else {
			for ( auto edge: set_inner_edges ) {
				contractions[edge] = nodes_contraction<ndims_type, size_type>() ;

				contractions.at( edge ).nodes2contract = current_edges2nodes.at( edge ) ;
				contractions.at( edge ).new_node = *current_edges2nodes.at( edge ).begin() ;
				std::set<ndims_type> edges_in_union ;
				for ( auto node: contractions.at( edge ).nodes2contract ) {
					for ( auto edge_in_union: current_nodes2edges.at( node ) ) {
						edges_in_union.insert( edge_in_union ) ;
					}
					if ( node < contractions.at( edge ).new_node ) { contractions.at( edge ).new_node = node ; }
				}
				for ( auto edge_in_union: edges_in_union ) {
					if ( set_outer_edges.count( edge_in_union ) ) {
						contractions.at( edge ).edges2keep.insert( edge_in_union ) ;
					}
					else {
						contractions.at( edge ).edges2contract.insert( edge_in_union ) ;
						for ( auto node: current_edges2nodes.at( edge_in_union ) ) {
							if ( ! contractions.at( edge ).nodes2contract.count( node ) ) {
								contractions.at( edge ).edges2contract.erase( edge_in_union ) ;
								contractions.at( edge ).edges2keep.insert( edge_in_union ) ;
								break ;
							}
						}
					}
				}
				contractions.at( edge ).removed_size = 1 ;
				for ( auto edge2contract: contractions.at( edge ).edges2contract ) {
					contractions.at( edge ).removed_size *= edges2sizes.at( edge2contract ) ;
				}
				contractions.at( edge ).cost = 1 ;
				for ( auto edge_in_union: edges_in_union ) {
					contractions.at( edge ).cost *= edges2sizes.at( edge_in_union ) ;
				}
			}
		}

		// find best contraction (i.e., greedy)
		ndims_type best_contraction = *set_inner_edges.begin() ;
		for ( auto edge: set_inner_edges ) {
			if ( contractions.at( edge ).removed_size > contractions.at( best_contraction ).removed_size ||
					( contractions.at( edge ).removed_size == contractions.at( best_contraction ).removed_size
							&& contractions.at( edge ).cost < contractions.at( best_contraction ).cost ) ) {
				best_contraction = edge ;
			}
		}

		// add to contraction_sequence
		contraction_sequence.push_back( contractions.at( best_contraction ) ) ;

		// assign nodes2edges, edges2nodes
		contraction_sequence.back().edges2nodes = current_edges2nodes ;
		contraction_sequence.back().nodes2edges = current_nodes2edges ;

		// update current_edges2nodes and set_inner_edges
		if ( contraction_sequence.back().nodes2contract.size() > 1 ) {
			for ( auto edge: contraction_sequence.back().edges2keep ) {
				for ( auto node: contraction_sequence.back().nodes2contract ) {
					current_edges2nodes.at( edge ).erase( node ) ;
				}
				current_edges2nodes.at( edge ).insert( contraction_sequence.back().new_node ) ;
			}
		}
		for ( auto edge: contraction_sequence.back().edges2contract ) {
			current_edges2nodes.erase( edge ) ;
			set_inner_edges.erase( edge ) ;
		}

		// update current_nodes2edges and current_nodes2shapes
		for ( auto node: contraction_sequence.back().nodes2contract ) {
			if ( node == contraction_sequence.back().new_node ) {
				current_nodes2edges.at( node ).clear() ;
			}
			else {
				current_nodes2edges.erase( node ) ;
			}
		}

		// assign new_nodes_edges
		if ( set_inner_edges.size() == 0 && current_nodes2edges.size() == 1 ) {
			contraction_sequence.back().new_node_edges = outer_edges ;
		}
		else {
			for ( auto edge: contraction_sequence.back().edges2keep ) {
				contraction_sequence.back().new_node_edges.push_back( edge ) ;
			}
		}
		for ( auto edge: contraction_sequence.back().new_node_edges ) {
			current_nodes2edges.at( contraction_sequence.back().new_node ).push_back( edge ) ;
		}

		// construct multi_indices
		std::vector<boost::shared_ptr<vector_wrapper<size_type>>> selection_arrays ;
		std::vector<size_type> shape_out_outer_index, shape_out_index ;
		std::map<ndims_type, ndims_type> outer_edges2entry_of_outer_index, inner_edges2entry_of_index ;
		std::map<ndims_type, boost::shared_ptr<dense_vector<size_type>>> nodes2index_sel ;

		j = 0 ;
		for ( auto edge: contraction_sequence.back().new_node_edges ) {
			shape_out_outer_index.push_back( edges2sizes.at( edge ) ) ;
			outer_edges2entry_of_outer_index[edge] = j ;
			++j ;
		}
		contraction_sequence.back().outer_index = boost::make_shared<shape_index<size_type>>( shape_out_outer_index ) ;
		contraction_sequence.back().new_node_shape = shape_out_outer_index ;
		for ( auto edge: contraction_sequence.back().new_node_edges ) {
			shape_out_index.push_back( 1 ) ;
			selection_arrays.push_back( boost::make_shared<vector_wrapper<size_type>>( (*contraction_sequence.back().outer_index)[dense_vector<size_type>( outer_edges2entry_of_outer_index[edge] )] ) ) ;
		}

		for ( auto node: contraction_sequence.back().nodes2contract ) {
			std::vector<size_type> dum ;
			for ( auto edge: contraction_sequence.back().nodes2edges.at( node ) ) {
				if ( contraction_sequence.back().edges2contract.count( edge ) ) {
					if ( ! inner_edges2entry_of_index.count( edge ) ) {
						shape_out_index.push_back( edges2sizes.at( edge ) ) ;
						selection_arrays.push_back( boost::make_shared<vector_wrapper<size_type>>( range< size_type >( 0, edges2sizes.at( edge ) ) ) ) ;
						inner_edges2entry_of_index[edge] = j ;
						++j ;
					}
					dum.push_back( inner_edges2entry_of_index[edge] ) ;
				}
				else {
					dum.push_back( outer_edges2entry_of_outer_index[edge] ) ;
				}
			}
			nodes2index_sel[node] = boost::make_shared<dense_vector<size_type>>( dum ) ;
		}

		contraction_sequence.back().index = boost::make_shared<block_selection_index<vector_wrapper<size_type>>>( shape_out_index, selection_arrays ) ;

		for ( auto nodes2index_sel_it: nodes2index_sel ) {
			contraction_sequence.back().nodes2indices[nodes2index_sel_it.first] =
					boost::make_shared<decltype( (*contraction_sequence.back().index)[*nodes2index_sel_it.second] )>
			        ( (*contraction_sequence.back().index)[*nodes2index_sel_it.second] ) ;
		}

	}

	return contraction_sequence ;
}

template <typename NE, typename EN, typename ES, typename IE, typename OE >
typename std::deque<nodes_contraction<typename IE::value_type, typename ES::mapped_type>> greedy_sequence_entry( NE nodes2edges, EN edges2nodes,
		ES edges2sizes_orig, IE inner_edges, OE outer_edges ) {

	typedef typename IE::value_type ndims_type ;
	typedef typename ES::mapped_type size_type ;

	ndims_type i, j, k ;

	std::deque<nodes_contraction<ndims_type, size_type>> contraction_sequence ;

	std::set<ndims_type> set_inner_edges( inner_edges.begin(), inner_edges.end() ) ;
	std::set<ndims_type> set_outer_edges( outer_edges.begin(), outer_edges.end() ) ;
	std::map<ndims_type, ndims_type> outer_edges2dim ;
	std::set<ndims_type> reduced_outer_edges ;
	boost::shared_ptr<dense_vector<size_type>> entry_index = boost::make_shared<dense_vector<size_type>>( 0, outer_edges.size() ) ;
	i = 0 ;
	for ( auto edge: outer_edges ) {
		outer_edges2dim[edge] = i ;
		++i ;
	}

	std::map<ndims_type, std::set<ndims_type>>         current_edges2nodes ;
	std::map<ndims_type, std::vector<ndims_type>>      current_nodes2edges ;
	std::map<ndims_type, size_type>                    edges2sizes ;
	for ( auto edges2nodes_it: edges2nodes ) {
		auto edge = edges2nodes_it.first ;
		for ( auto nodes_it: edges2nodes_it.second ) {
			current_edges2nodes[edge].insert( std::get<0>( nodes_it ) ) ;
		}
	}
	for ( auto nodes2edges_it: nodes2edges ) {
		auto node = nodes2edges_it.first ;
		for ( auto edge: nodes2edges_it.second ) {
			current_nodes2edges[node].push_back( edge ) ;
		}
	}
	for ( auto edges2sizes_it: edges2sizes_orig ) {
		auto edge = edges2sizes_it.first ;
		if ( set_outer_edges.count( edge ) ) {
			edges2sizes[edge] = 1 ;
		}
		else {
			edges2sizes[edge] = edges2sizes_it.second ;
		}
	}

	while ( set_inner_edges.size() > 0 || current_nodes2edges.size() > 1 ) {
		std::map<ndims_type, nodes_contraction<ndims_type, size_type>> contractions ;
		if ( set_inner_edges.size() == 0 ) {
			// store in local variable
			contractions[0].new_node = current_nodes2edges.begin()->first ;
			for ( auto node2edges: current_nodes2edges ) {
				contractions.at( 0 ).nodes2contract.insert( node2edges.first ) ;
				if ( node2edges.first < contractions.at( 0 ).new_node ) { contractions.at( 0 ).new_node = node2edges.first ; }
			}
			for ( auto node: contractions.at( 0 ).nodes2contract ) {
				for ( auto edge: current_nodes2edges.at( node ) ) {
					contractions.at( 0 ).edges2keep.insert( edge ) ;
				}
			}
		}
		else {
			for ( auto edge: set_inner_edges ) {
				contractions[edge] = nodes_contraction<ndims_type, size_type>() ;

				contractions.at( edge ).nodes2contract = current_edges2nodes.at( edge ) ;
				contractions.at( edge ).new_node = *current_edges2nodes.at( edge ).begin() ;
				std::set<ndims_type> edges_in_union ;
				for ( auto node: contractions.at( edge ).nodes2contract ) {
					for ( auto edge_in_union: current_nodes2edges.at( node ) ) {
						edges_in_union.insert( edge_in_union ) ;
					}
					if ( node < contractions.at( edge ).new_node ) { contractions.at( edge ).new_node = node ; }
				}
				for ( auto edge_in_union: edges_in_union ) {
					if ( set_outer_edges.count( edge_in_union ) ) {
						contractions.at( edge ).edges2keep.insert( edge_in_union ) ;
					}
					else {
						contractions.at( edge ).edges2contract.insert( edge_in_union ) ;
						for ( auto node: current_edges2nodes.at( edge_in_union ) ) {
							if ( ! contractions.at( edge ).nodes2contract.count( node ) ) {
								contractions.at( edge ).edges2contract.erase( edge_in_union ) ;
								contractions.at( edge ).edges2keep.insert( edge_in_union ) ;
								break ;
							}
						}
					}
				}
				contractions.at( edge ).removed_size = 1 ;
				for ( auto edge2contract: contractions.at( edge ).edges2contract ) {
					contractions.at( edge ).removed_size *= edges2sizes.at( edge2contract ) ;
				}
				contractions.at( edge ).cost = 1 ;
				for ( auto edge_in_union: edges_in_union ) {
					contractions.at( edge ).cost *= edges2sizes.at( edge_in_union ) ;
				}
			}
		}

		// find best contraction (i.e., greedy)
		ndims_type best_contraction = *set_inner_edges.begin() ;
		for ( auto edge: set_inner_edges ) {
			if ( contractions.at( edge ).removed_size > contractions.at( best_contraction ).removed_size ||
					( contractions.at( edge ).removed_size == contractions.at( best_contraction ).removed_size
							&& contractions.at( edge ).cost < contractions.at( best_contraction ).cost ) ) {
				best_contraction = edge ;
			}
		}

		// add to contraction_sequence
		contraction_sequence.push_back( contractions.at( best_contraction ) ) ;

		// assign nodes2edges, edges2nodes
		contraction_sequence.back().edges2nodes = current_edges2nodes ;
		contraction_sequence.back().nodes2edges = current_nodes2edges ;

		// update current_edges2nodes and set_inner_edges
		if ( contraction_sequence.back().nodes2contract.size() > 1 ) {
			for ( auto edge: contraction_sequence.back().edges2keep ) {
				for ( auto node: contraction_sequence.back().nodes2contract ) {
					current_edges2nodes.at( edge ).erase( node ) ;
				}
				current_edges2nodes.at( edge ).insert( contraction_sequence.back().new_node ) ;
			}
		}
		for ( auto edge: contraction_sequence.back().edges2contract ) {
			current_edges2nodes.erase( edge ) ;
			set_inner_edges.erase( edge ) ;
		}

		// update current_nodes2edges and current_nodes2shapes
		for ( auto node: contraction_sequence.back().nodes2contract ) {
			if ( node == contraction_sequence.back().new_node ) {
				current_nodes2edges.at( node ).clear() ;
			}
			else {
				current_nodes2edges.erase( node ) ;
			}
		}

		// assign new_nodes_edges
		if ( set_inner_edges.size() == 0 && current_nodes2edges.size() == 1 ) {
			contraction_sequence.back().new_node_edges = outer_edges ;
		}
		else {
			for ( auto edge: contraction_sequence.back().edges2keep ) {
				contraction_sequence.back().new_node_edges.push_back( edge ) ;
			}
		}
		for ( auto edge: contraction_sequence.back().new_node_edges ) {
			current_nodes2edges.at( contraction_sequence.back().new_node ).push_back( edge ) ;
		}

		// assign entry_index
		contraction_sequence.back().entry_index = entry_index ;

		// construct multi_indices
		std::vector<boost::shared_ptr<vector_wrapper<size_type>>> selection_arrays ;
		std::vector<size_type> shape_out_outer_index, shape_out_index ;
		std::map<ndims_type, ndims_type> outer_edges2entry_of_outer_index, inner_edges2entry_of_index ;
		std::map<ndims_type, boost::shared_ptr<dense_vector<size_type>>> nodes2index_sel ;

		j = 0 ;
		for ( auto edge: contraction_sequence.back().new_node_edges ) {
			shape_out_outer_index.push_back( edges2sizes.at( edge ) ) ;
			outer_edges2entry_of_outer_index[edge] = j ;
			++j ;
		}
		contraction_sequence.back().outer_index = boost::make_shared<shape_index<size_type>>( shape_out_outer_index ) ;
		contraction_sequence.back().new_node_shape = shape_out_outer_index ;
		for ( auto edge: contraction_sequence.back().new_node_edges ) {
			shape_out_index.push_back( 1 ) ;
			if ( set_outer_edges.count( edge ) && ! reduced_outer_edges.count( edge ) ) {
				selection_arrays.push_back( boost::make_shared<vector_wrapper<size_type>>
						( (*contraction_sequence.back().entry_index)[dense_vector<size_type>( outer_edges2dim[edge] )] ) ) ;
				reduced_outer_edges.insert( edge ) ;
			}
			else {
				selection_arrays.push_back( boost::make_shared<vector_wrapper<size_type>>
						( (*contraction_sequence.back().outer_index)[dense_vector<size_type>( outer_edges2entry_of_outer_index[edge] )] ) ) ;
			}
		}

		for ( auto node: contraction_sequence.back().nodes2contract ) {
			std::vector<size_type> dum ;
			for ( auto edge: contraction_sequence.back().nodes2edges.at( node ) ) {
				if ( contraction_sequence.back().edges2contract.count( edge ) ) {
					if ( ! inner_edges2entry_of_index.count( edge ) ) {
						shape_out_index.push_back( edges2sizes.at( edge ) ) ;
						selection_arrays.push_back( boost::make_shared<vector_wrapper<size_type>>( range< size_type >( 0, edges2sizes.at( edge ) ) ) ) ;
						inner_edges2entry_of_index[edge] = j ;
						++j ;
					}
					dum.push_back( inner_edges2entry_of_index[edge] ) ;
				}
				else {
					dum.push_back( outer_edges2entry_of_outer_index[edge] ) ;
				}
			}
			nodes2index_sel[node] = boost::make_shared<dense_vector<size_type>>( dum ) ;
		}

		contraction_sequence.back().index = boost::make_shared<block_selection_index<vector_wrapper<size_type>>>( shape_out_index, selection_arrays ) ;

		for ( auto nodes2index_sel_it: nodes2index_sel ) {
			contraction_sequence.back().nodes2indices[nodes2index_sel_it.first] =
					boost::make_shared<decltype( (*contraction_sequence.back().index)[*nodes2index_sel_it.second] )>
			        ( (*contraction_sequence.back().index)[*nodes2index_sel_it.second] ) ;
		}

	}

	return contraction_sequence ;
}

} // namespace glas3

#endif

