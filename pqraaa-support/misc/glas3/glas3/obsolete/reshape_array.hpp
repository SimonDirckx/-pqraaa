//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_reshape_array_hpp
#define glas3_array_dense_array_type_reshape_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/concept/array.hpp>

#include <glas3/array/algorithm/indexing.hpp>

#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

namespace glas3 {

template <typename S>
class reshape_array {
    static_assert( is<Array, S>::value, "S should be an Array" ) ;

    public:
      typedef std::ptrdiff_t                     value_type;
      typedef S                                  shape_type ;
      typedef typename shape_type::value_type    size_type ;
      typedef typename shape_type::size_type     ndims_type ;

    public:
      reshape_array () // empty reshape_array
      : size_( 1 )
      {}

      explicit reshape_array ( shape_type shape )
      : shape_( boost::make_shared< shape_type >( std::move( shape ) ) )
      , size_( 1 )
      {
  	    for ( ndims_type k = 0; k < shape_->size(); ++k ) {
  		    size_ *= (*shape_)[k] ;
  	    }
      }

    public:
        // Copy constructor -> do not allow
        reshape_array ( reshape_array const& that ) = delete ;

        // Move constructor
        reshape_array ( reshape_array && that ) = default;

    public:
        // Copy assignment -> do not allow
        reshape_array& operator= ( reshape_array const& that ) = delete ; // not assignable

        // Move assignment
        reshape_array& operator= ( reshape_array&& that) = default;

    private:
        // Constructor used in shallow_copy
        explicit reshape_array ( boost::shared_ptr< shape_type > shape, size_type size )
        : shape_( shape )
        , size_( size )
        {}

    private:
        boost::shared_ptr< shape_type >    shape_ ;
        size_type                          size_ ;

    public:
        reshape_array shallow_copy () const {
           	return reshape_array( shape_, size_ );
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
            assert( i >= 0 && i < size_ ) ; // instead of this check apply mod operator?
            return i ;
        }

    	template <typename X>
    	typename std::enable_if< is<Vector, X>::value, value_type >::type
    	operator() ( X const& j ) const {
    		assert( j.size() == shape_->size() ) ;
    		return multi_index2index_initializer_list( *shape_, j );
    	}

  } ;


template <typename S>
struct concept< reshape_array<S>, typename std::enable_if< is<Array, S>::value >::type >
: DenseArray
{};

} // namespace glas3


#endif
