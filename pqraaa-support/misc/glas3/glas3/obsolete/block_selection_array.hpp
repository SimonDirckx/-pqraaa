//  (C) Copyright Sam Corveleyn 2015.
//  Use, modification and distribution are subject to the
//  GLAS Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt)

#ifndef glas3_array_dense_array_type_slice_selection_array_hpp
#define glas3_array_dense_array_type_slice_selection_array_hpp

#include <glas3/concept/is.hpp>
#include <glas3/concept/concept.hpp>
#include <glas3/array/dense_array/concept/dense_array.hpp>

#include <glas3/array/dense_array/container/dense_vector.hpp>

#include <initializer_list>
#include <vector>
#include <type_traits>
#include <cassert>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

#include <iostream>

namespace glas3 {

// WARNING: slice_selection_array makes shallow_copy of selection arrays. Changing the passed selection arrays afterwards results, however, in undefined behavior,
// except if you reset the internal index memory by executing .resetIndexMemory()

template <typename X>
class slice_selection_array {
	static_assert( is<Array, X>::value, "V should be integral" ) ;

    public:
        typedef typename X::value_type            value_type ;
        typedef dense_vector<std::ptrdiff_t>      shape_type ;
        typedef typename shape_type::value_type   size_type ;
        typedef typename shape_type::size_type    ndims_type ;

    public:
        explicit slice_selection_array ()
        : selection_arrays_( )
        , shape_orig_( boost::make_shared< shape_type >( 1 ) )
        , shape_( boost::make_shared< shape_type >( 1 ) )
        , dim_factor_orig_( boost::make_shared< shape_type >( 1 ) )
        , dim_factor_( boost::make_shared< shape_type >( 1 ) )
        , j_mem_ ( )
        , i_in_ ( 0 )
        , i_out_ ( 0 )
        , size_( 1 )
        , ndof_( 1 )
    	{}

        // Copy constructor from std::initializer_list< selection_array_type > -> shallow copy
        template < typename S >
        slice_selection_array ( S const& shape_orig, std::vector< X > const& selection_array_list, typename std::enable_if<is<Array, S>::value>::type* = 0 )
        : selection_arrays_( selection_array_list.size() )
        , shape_orig_( boost::make_shared< shape_type >( 0, shape_orig.size() ) )
        , shape_( boost::make_shared< shape_type >( 0, shape_orig.size() ) )
        , dim_factor_orig_( boost::make_shared< shape_type >( 1, shape_orig.size() ) )
        , dim_factor_( boost::make_shared< shape_type >( 1, shape_orig.size() ) )
        , j_mem_ ( boost::make_shared< shape_type >( 0, shape_orig.size() ) )
        , i_in_ ( boost::make_shared< size_type >( 0 ) )
        , i_out_ ( boost::make_shared< size_type >( 0 ) )
        , size_( 1 )
        , ndof_( 0 )
        {
        	assert( shape_orig.size() == selection_array_list.size() ) ;
        	for ( ndims_type k = 0; k < shape_orig.size(); ++k ){
        		selection_arrays_[k] = boost::make_shared< X >( std::move( selection_array_list[k].shallow_copy() ) ) ;
        		(*shape_orig_)[k] = shape_orig[k] ;
        		(*shape_)[k] = selection_arrays_[k]->size() ;
        		size_ *= (*shape_)[k] ;
        		ndof_ += selection_arrays_[k]->ndof() ;
        		if ( k > 0 ) {
        			(*dim_factor_orig_)[k] = (*dim_factor_orig_)[k - 1] * (*shape_orig_)[k - 1] ;
        			(*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ;
        		}
        		*i_out_ += (*(selection_arrays_[k]))[0] * (*dim_factor_orig_)[k] ;
        	}
        }

        template < typename S >
        slice_selection_array ( S const& shape_orig, std::vector<boost::shared_ptr<X>> const& selection_array_list, typename std::enable_if<is<Array, S>::value>::type* = 0 )
        : selection_arrays_( selection_array_list )
        , shape_orig_( boost::make_shared< shape_type >( 0, shape_orig.size() ) )
        , shape_( boost::make_shared< shape_type >( 0, shape_orig.size() ) )
        , dim_factor_orig_( boost::make_shared< shape_type >( 1, shape_orig.size() ) )
        , dim_factor_( boost::make_shared< shape_type >( 1, shape_orig.size() ) )
        , j_mem_ ( boost::make_shared< shape_type >( 0, shape_orig.size() ) )
        , i_in_ ( boost::make_shared< size_type >( 0 ) )
        , i_out_ ( boost::make_shared< size_type >( 0 ) )
        , size_( 1 )
        , ndof_( 0 )
        {
        	assert( shape_orig.size() == selection_array_list.size() ) ;
        	for ( ndims_type k = 0; k < shape_orig.size(); ++k ){
        		//selection_arrays_[k] = boost::make_shared< X >( std::move( (*(selection_array_list[k])).shallow_copy() ) ) ;
        		(*shape_orig_)[k] = shape_orig[k] ;
        		(*shape_)[k] = selection_arrays_[k]->size() ;
        		size_ *= (*shape_)[k] ;
        		ndof_ += selection_arrays_[k]->ndof() ;
        		if ( k > 0 ) {
        			(*dim_factor_orig_)[k] = (*dim_factor_orig_)[k - 1] * (*shape_orig_)[k - 1] ;
        			(*dim_factor_)[k] = (*dim_factor_)[k - 1] * (*shape_)[k - 1] ;
        		}
        		*i_out_ += (*(selection_arrays_[k]))[0] * (*dim_factor_orig_)[k] ;
        	}
        }

    public:
	    // Copy constructor -> do not allow
        slice_selection_array ( slice_selection_array const& that ) = delete ;

	    // Move constructor
        slice_selection_array ( slice_selection_array && ) = default;

    public:
        // Copy assignment -> do not allow
        slice_selection_array& operator= ( slice_selection_array const& that ) = delete ; // not assignable

        // Move assignment
        slice_selection_array& operator= ( slice_selection_array&& that ) = default;

    private:
        // Constructor used in shallow_copy
        explicit slice_selection_array ( std::vector<boost::shared_ptr<X>> selection_arrays, boost::shared_ptr<shape_type> shape_orig, boost::shared_ptr<shape_type> shape,
        		boost::shared_ptr<shape_type> dim_factor_orig, boost::shared_ptr<shape_type> dim_factor, boost::shared_ptr<shape_type> j_mem, boost::shared_ptr<size_type> i_in,
        		boost::shared_ptr<size_type> i_out, size_type size, size_type ndof )
        : selection_arrays_( selection_arrays )
        , shape_orig_( shape_orig )
        , shape_( shape )
        , dim_factor_orig_( dim_factor_orig )
        , dim_factor_( dim_factor )
        , j_mem_ ( j_mem )
        , i_in_ ( i_in )
        , i_out_ ( i_out )
        , size_( size )
        , ndof_( ndof )
        {}

    private:
        std::vector<boost::shared_ptr<X>>                     selection_arrays_ ;
        boost::shared_ptr<shape_type>                         shape_orig_ ;
        boost::shared_ptr<shape_type>                         shape_ ;
        boost::shared_ptr<shape_type>                         dim_factor_orig_ ;
        boost::shared_ptr<shape_type>                         dim_factor_ ;
        boost::shared_ptr<shape_type>                         j_mem_ ;
        boost::shared_ptr<size_type>                          i_in_ ;
        boost::shared_ptr<size_type>                          i_out_ ;
        size_type                                             size_ ;
        size_type                                             ndof_ ;

    public:
    	slice_selection_array shallow_copy () const {
           	return slice_selection_array( selection_arrays_, shape_orig_, shape_, dim_factor_orig_, dim_factor_, j_mem_, i_in_, i_out_, size_, ndof_ );
        }

    public:
        shape_type const& shape () const {
            return *shape_ ;
        }

        size_type size () const {
            return size_ ;
        }

        size_type ndof () const {
            return ndof_ ;
        }

    public:
        void resetIndexMemory() const {
        	*i_out_ = 0 ;
        	for ( ndims_type k = 0; k < selection_arrays_.size(); ++k ){
        		*i_out_ += (*(selection_arrays_[k]))[(*j_mem_)[k]] * (*dim_factor_orig_)[k] ;
        	}
        }

    public:
        value_type
        operator[] ( size_type const& i ) const {
        	assert( i >= 0 && i < size_ ) ;
        	index2multi_index2index_selection_arrays_memory( *dim_factor_orig_, *shape_orig_, *dim_factor_, *shape_, selection_arrays_, *j_mem_, i, *i_in_, *i_out_ ) ;
        	return *i_out_ ;
        }

    	template <typename S>
    	typename std::enable_if< is<Vector, S>::value, value_type >::type
    	operator() ( S const& j ) const {
    		assert( j.size() == shape_->size() ) ;
    		return multi_index2index_selection_arrays_initializer_list( *shape_orig_, *shape_, selection_arrays_, j );
    	}

} ;

template <typename X>
struct concept< slice_selection_array<X>, typename std::enable_if< is<Array, X>::value >::type >
: DenseArray
{};

} ; // namespace glas3

#endif
