//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_permute_array_hpp
#define glas3_array_dense_array_type_permute_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <glas3/array/algorithm/indexing.hpp>

#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

namespace glas3 {

template <typename S, typename D>
class permute_array {
	static_assert( is<Array, S>::value, "S should be an Array" ) ;
	static_assert( is<Array, D>::value, "D should be an Array" ) ;

    public:
      typedef typename S::value_type                  value_type ;
      typedef dense_vector<value_type>                shape_type ;
      typedef typename shape_type::value_type         size_type ;
      typedef typename shape_type::size_type          ndims_type ;

    public:
        permute_array () // empty permute_array
        : i_in_ ( boost::make_shared< size_type >( 0 ) )
        , i_out_ ( boost::make_shared< size_type >( 0 ) )
        , size_ ( 1 )
        {}

        explicit permute_array ( S const& shape_orig, D dim_order )
        : shape_( boost::make_shared<shape_type> ( 0, shape_orig.size() ) )
        , dim_order_( boost::make_shared<D> ( std::move( dim_order ) ) )
        , dim_factor_( boost::make_shared<shape_type> ( 1, shape_orig.size() ) )
        , dim_factor_asifnotpermuted_( boost::make_shared<shape_type> ( 1, shape_orig.size() ) )
        , j_mem_ ( boost::make_shared< shape_type >( 0, shape_orig.size() ) )
        , i_in_ ( boost::make_shared< size_type >( 0 ) )
        , i_out_ ( boost::make_shared< size_type >( 0 ) )
        , size_ ( 1 )
        {
    	    assert( shape_orig.size() == dim_order_->size() ) ;

    	    shape_type dim_factor_orig( 1, shape_orig.size() ) ;

    	    for ( ndims_type k = 0; k < shape_->size(); ++k ) {
    		    size_ *= shape_orig[k] ;
    		    (*shape_)[k] = shape_orig[(*dim_order_)[k]] ;
    		    if ( k > 0 ) {
    			    dim_factor_orig[k] = dim_factor_orig[k - 1] * shape_orig[k - 1];
    			    (*dim_factor_asifnotpermuted_)[k] = (*dim_factor_asifnotpermuted_)[k - 1] * (*shape_)[k - 1];
    		    }
    	    }

    	    for ( ndims_type k = 0; k < shape_->size(); ++k ) {
    		    (*dim_factor_)[k] = dim_factor_orig[(*dim_order_)[k]] ;
    	    }
        }

    public:
        // Copy constructor -> do not allow
        permute_array ( permute_array const& that ) = delete ;

        // Move constructor
        permute_array ( permute_array && that ) = default;

    public:
        // Copy assignment -> do not allow
        permute_array& operator= ( permute_array const& that ) = delete ; // not assignable

        // Move assignment
        permute_array& operator= ( permute_array&& that) = default;

    private:
        // Constructor used in shallow_copy
        explicit permute_array ( boost::shared_ptr<shape_type> shape, boost::shared_ptr<D> dim_order, boost::shared_ptr<shape_type> dim_factor,
        		boost::shared_ptr<shape_type> dim_factor_asifnotpermuted, boost::shared_ptr< shape_type > j_mem,
        		boost::shared_ptr< size_type > i_in, boost::shared_ptr< size_type > i_out, size_type size )
        : shape_( shape )
        , dim_order_ ( dim_order )
        , dim_factor_ ( dim_factor )
        , dim_factor_asifnotpermuted_ ( dim_factor_asifnotpermuted )
        , j_mem_ ( j_mem )
        , i_in_ ( i_in )
        , i_out_ ( i_out )
        , size_( size )
        {}

    private:
        boost::shared_ptr< shape_type >     shape_ ;
        boost::shared_ptr< D >              dim_order_ ;
        boost::shared_ptr< shape_type >     dim_factor_ ;
        boost::shared_ptr< shape_type >     dim_factor_asifnotpermuted_ ;
        boost::shared_ptr< shape_type >     j_mem_ ;
        boost::shared_ptr< size_type >      i_in_ ;
        boost::shared_ptr< size_type >      i_out_ ;
        size_type                           size_ ;

    public:
        permute_array shallow_copy () const {
           	return permute_array( shape_, dim_order_, dim_factor_, dim_factor_asifnotpermuted_, j_mem_, i_in_, i_out_, size_ );
        }

    public:
        shape_type const& shape () const {
            return *shape_ ;
        }

        size_type size () const {
            return size_ ;
        }

        size_type ndof () const {
            return 0 ;
        }

    public:
        value_type
        operator[] ( size_type const& i ) const {
            assert( i >= 0 && i < size_ ) ;
            index2multi_index2index_memory( *dim_factor_, *dim_factor_asifnotpermuted_, *shape_, *j_mem_, i, *i_in_, *i_out_) ;
            return *i_out_ ;
        }

    	template <typename X>
    	typename std::enable_if< is<Vector, X>::value, value_type >::type
    	operator() ( X const& j ) const {
    		assert( j.size() == shape_->size() ) ;
    		return multi_index2index_dim_factor_initializer_list( *shape_, *dim_factor_, j );
    	}

  } ;


  template <typename S, typename D>
  struct concept< permute_array<S, D>, typename std::enable_if< is<Array, S>::value && is<Array, D>::value >::type >
  : DenseArray
  {};

} // namespace glas3


#endif
