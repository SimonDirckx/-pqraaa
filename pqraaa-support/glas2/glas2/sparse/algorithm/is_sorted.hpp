#ifndef glas2_sparse_algorithm_is_sorted_hpp
#define glas2_sparse_algorithm_is_sorted_hpp

namespace glas2 {

  template <typename C>
  bool is_sorted( C const& coo) {
    bool sorted = true ;

    if (coo.num_nz()==0) return true;

    typename C::size_type index = 1 ;
    typename C::size_type previous_row = coo.row(0);
    typename C::size_type previous_column = coo.column(0);

    while (sorted && index<coo.num_nz()) {
      if (coo.row(index) < previous_row) sorted=false;
      if (coo.row(index) == previous_row && coo.column(index)<=previous_column) sorted=false;

      previous_row = coo.row(index) ;
      previous_column = coo.column(index) ;
      ++index ;
    }

    return sorted ;
  }

}

#endif
