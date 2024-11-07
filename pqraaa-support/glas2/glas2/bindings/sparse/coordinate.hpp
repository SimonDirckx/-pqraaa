//  (C) Copyright Karl Meerbergen, 2011.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_bindings_sparse_coordinate_hpp
#define glas2_bindings_sparse_coordinate_hpp

#include <glas2/concept/is.hpp>
#include <glas2/sparse/concept/coordinate_sparse_matrix.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/tag.hpp>
#include <type_traits>
#include <tuple>
#include <utility>


namespace boost { namespace numeric { namespace bindings { namespace detail {

  template< typename M, typename Id >
  struct adaptor< M, Id, typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, M >::value >::type >
  {
    typedef typename std::conditional< std::is_const< Id >::value
                                     , typename M::value_type const
                                     , typename M::value_type
                                     >::type value_type;
    typedef typename std::conditional< std::is_const< Id >::value
                                     , typename M::index_type const
                                     , typename M::index_type
                                     >::type index_type;

    typedef std::ptrdiff_t size_type ;
    //typedef typename M::size_type size_type ;

    typedef std::tuple<
        std::pair< tag::value_type, value_type >,
        std::pair< tag::index_type, index_type >,
        std::pair< tag::entity, tag::matrix >,
        std::pair< tag::size_type<1>, typename M::size_type >,
        std::pair< tag::size_type<2>, typename M::size_type >,
        std::pair< tag::data_structure, tag::coordinate_sparse >,
        std::pair< tag::matrix_type, tag::general >,  // At the moment only general nonsymmetric matrices are supported.
        std::pair< tag::data_order, tag::row_major >, // Assumption although this is not guaranteed
        std::pair< tag::index_base, std::integral_constant< int, Id::index_base > >
    > property_map ;

    static size_type size1( const M& id ) {
        return id.num_rows();
    }

    static size_type size2( const M& id ) {
        return id.num_columns();
    }

    static value_type* begin_value( Id& id ) {
        return id.data().ptr() ;
    }

    static value_type* end_value( Id& id ) {
        return id.data().ptr() + id.num_nz() ;
    }

/*    static std::ptrdiff_t stride1( const Id& id ) {
        return glas::stride_1( id ) ;
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return glas::stride_2( id ) ;
    }
*/
    static index_type* begin_index_major( Id& id ) {
        return id.row_indices().ptr() ;
    }

    static index_type* end_index_major( Id& id ) {
        return id.row_indices().ptr() + id.num_nz() ;
    }

    static index_type* begin_index_minor( Id& id ) {
        return id.column_indices().ptr() ;
    }

    static index_type* end_index_minor( Id& id ) {
        return id.column_indices().ptr() + id.num_nz() ;
    }

  } ;



} } } } // boost::numeric::bindings::detail

#endif
