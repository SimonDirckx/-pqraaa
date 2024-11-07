//  (C) Copyright Sam Corveleyn 2014.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_tile_array_hpp
#define glas3_array_dense_array_type_tile_array_hpp

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

template <typename S_Mi, typename S_Ma>
class tile_array {
	static_assert( is<Array, S_Mi>::value, "S_Mi should be an Array" ) ;
	static_assert( is<Array, S_Ma>::value, "S_Ma should be an Array" ) ;

    public:
      typedef std::ptrdiff_t                          value_type ;
      typedef dense_vector<value_type>                shape_type ;
      typedef typename shape_type::value_type         size_type ;
      typedef typename shape_type::size_type          ndims_type ;

    public:
      tile_array () // empty permute_array
        : i_in_ ( boost::make_shared< size_type >( 0 ) )
        , i_out_ ( boost::make_shared< size_type >( 0 ) )
        , size_ ( 1 )
        {}

        explicit tile_array ( S_Mi micro_shape, S_Ma macro_shape )
        : shape_( boost::make_shared<shape_type> ( 0, micro_shape.size() ) )
        , micro_shape_( boost::make_shared<S_Mi> ( std::move( micro_shape ) ) )
        , macro_shape_( boost::make_shared<S_Ma> ( std::move( macro_shape ) ) )
        , dim_factor_( boost::make_shared<shape_type> ( 1, micro_shape_->size() ) )
        , micro_dim_factor_( boost::make_shared<shape_type> ( 1, micro_shape_->size() ) )
        , j_mem_ ( boost::make_shared< shape_type >( 0, micro_shape_->size() ) )
        , i_in_ ( boost::make_shared< size_type >( 0 ) )
        , i_out_ ( boost::make_shared< size_type >( 0 ) )
        , size_ ( 1 )
        {
    	    assert( micro_shape_->size() == macro_shape_->size() ) ;

    	    for ( ndims_type k = 0; k < shape_->size(); ++k ) {
    	    	(*shape_)[k] = (*micro_shape_)[k] * (*macro_shape_)[k] ;
    	    	size_ *= (*shape_)[k] ;
    		    if ( k > 0 ) {
    			    (*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1];
    			    (*micro_dim_factor_)[k] = (*micro_dim_factor_)[k - 1] * (*micro_shape_)[k - 1];
    		    }
    	    }
        }

    public:
        // Copy constructor -> do not allow
        tile_array ( tile_array const& that ) = delete ;

        // Move constructor
        tile_array ( tile_array && that ) = default;

    public:
        // Copy assignment -> do not allow
        tile_array& operator= ( tile_array const& that ) = delete ; // not assignable

        // Move assignment
        tile_array& operator= ( tile_array&& that ) = default;

    private:
        // Constructor used in shallow_copy
        explicit tile_array ( boost::shared_ptr<shape_type> shape, boost::shared_ptr< S_Mi > micro_shape, boost::shared_ptr< S_Ma > macro_shape,
        		boost::shared_ptr<shape_type> dim_factor, boost::shared_ptr<shape_type> micro_dim_factor, boost::shared_ptr< shape_type > j_mem,
        		boost::shared_ptr< size_type > i_in, boost::shared_ptr< size_type > i_out, size_type size )
        : shape_( shape )
        , micro_shape_ ( micro_shape )
        , macro_shape_ ( macro_shape )
        , dim_factor_ ( dim_factor )
        , micro_dim_factor_ ( micro_dim_factor )
        , j_mem_ ( j_mem )
        , i_in_ ( i_in )
        , i_out_ ( i_out )
        , size_( size )
        {}

    private:
        boost::shared_ptr< shape_type >     shape_ ;
        boost::shared_ptr< S_Mi >           micro_shape_ ;
        boost::shared_ptr< S_Ma >           macro_shape_ ;
        boost::shared_ptr< shape_type >     dim_factor_ ;
        boost::shared_ptr< shape_type >     micro_dim_factor_ ;
        boost::shared_ptr< shape_type >     j_mem_ ;
        boost::shared_ptr< size_type >      i_in_ ;
        boost::shared_ptr< size_type >      i_out_ ;
        size_type                           size_ ;

    public:
        tile_array shallow_copy () const {
           	return tile_array( shape_, micro_shape_, macro_shape_, dim_factor_, micro_dim_factor_, j_mem_, i_in_, i_out_, size_ );
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
            index2multi_index2index_mod_shape_memory( *micro_dim_factor_, *micro_shape_, *dim_factor_, *shape_, *j_mem_, i, *i_in_, *i_out_ ) ;
            return *i_out_ ;
        }

    	template <typename X>
    	typename std::enable_if< is<Vector, X>::value, value_type >::type
    	operator() ( X const& j ) const {
    		assert( j.size() == shape_->size() ) ;
    		return multi_index2index_mod_shape_initializer_list( *micro_shape_, *shape_, j );
    	}

  } ;


  template <typename S_Mi, typename S_Ma>
  struct concept< tile_array<S_Mi, S_Ma>, typename std::enable_if< is<Array, S_Mi>::value && is<Array, S_Ma>::value >::type >
  : DenseArray
  {};

} // namespace glas3


#endif
