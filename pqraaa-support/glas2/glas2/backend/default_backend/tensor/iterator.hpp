//  (C) Copyright Karl Meerbergen 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_backend_default_backend_tensor_iterator_hpp
#define glas2_backend_default_backend_tensor_iterator_hpp

#include <glas2/vector/container/shared_vector.hpp>
#include <glas2/tensor/concept/dense_tensor.hpp>
#include <glas2/vector/concept/dense_vector.hpp>
#include <glas2/vector/algorithm/fill.hpp>
#include <glas2/concept/is.hpp>
#include <type_traits>
#include <cassert>

namespace glas2 { namespace tensor_detail {


  template <typename S>
  class iterator {
    public:
      typedef S                        shape_type ;
      typedef typename S::size_type    size_type ;
      typedef shared_vector<size_type> index_type ;

    public:
      iterator( shape_type shape )
      : shape_( shape )
      , index_( shape_.size()+1 )
      {
        fill( index_, 0 ) ;
      }

    public:
      bool is_end() {
        return index_( shape_.size() )>0 ;
      }

      typename vector_selection< typename index_type::base_type, range>::result_type operator*() const { return index_( range(0,shape_.size() ) ) ; }

      void operator++() {
        ++index_(0) ;
        int j = 0 ;
        bool cont = index_(j)>=shape_(j) ;
        while (cont) {
          index_(j) = 0 ;
          ++j ;
          ++index_(j) ;
          if (j>=shape_.size()) cont = false ;
          else cont = index_(j)>=shape_(j) ;
        }
      }

    private:
      shape_type const&  shape_ ;
      mutable index_type index_ ;
  } ;

} }  // namespace glas2::tensor_detail

#endif
