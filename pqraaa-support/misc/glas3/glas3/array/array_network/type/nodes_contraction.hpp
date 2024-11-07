//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_array_network_type_nodes_contraction_hpp
#define glas3_array_array_network_type_nodes_contraction_hpp

#include <glas3/array/type/vector_wrapper.hpp>
#include <glas3/array/dense_array/type/shape_index.hpp>
#include <glas3/array/dense_array/type/block_selection_index.hpp>

#include <set>
#include <vector>
#include <type_traits>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

namespace glas3 {

  template <typename N, typename S, typename EnableIf=void>
  struct nodes_contraction
  {} ;

  template <typename N, typename S>
  struct nodes_contraction<N, S, typename std::enable_if< std::is_integral<N>::value >::type >
  {
	  typedef N    ndims_type ;
	  typedef S    size_type ;

	  std::map<ndims_type, std::set<ndims_type>>                                                   edges2nodes ;
	  std::map<ndims_type, std::vector<ndims_type>>                                                nodes2edges ;
	  std::set<ndims_type>                                                                         nodes2contract ;
	  std::set<ndims_type>                                                                         edges2contract ;
	  std::set<ndims_type>                                                                         edges2keep ;
      size_type                                                                                    removed_size ;
	  size_type                                                                                    cost ;
	  boost::shared_ptr<block_selection_index<vector_wrapper<size_type>>>                          index ;
	  boost::shared_ptr<shape_index<size_type>>                                                    outer_index ;
	  std::map<ndims_type, boost::shared_ptr<decltype( (*index)[dense_vector<size_type>()] )>>     nodes2indices ;
	  boost::shared_ptr<dense_vector<size_type>>                                                   entry_index ;
	  ndims_type                                                                                   new_node ;
	  std::vector<size_type>                                                                       new_node_shape ;
	  std::vector<size_type>                                                                       new_node_edges ;
  } ;

  // definition of comparison function to circumvent incomprehensible bug
  template <typename ndims_type, typename size_type>
  bool operator==( nodes_contraction<ndims_type, size_type> const& a1, nodes_contraction<ndims_type, size_type> const& a2 ) {
      return false ;
  } ;


} // namespace glas3


#endif
