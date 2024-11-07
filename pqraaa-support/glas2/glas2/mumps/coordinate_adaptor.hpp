//  (C) Copyright Karl Meerbergen, 2011.
//  Use, modification and distribution are subject to the 
//  GLAS Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef glas2_mumps_coordinate_adaptor_hpp
#define glas2_mumps_coordinate_adaptor_hpp

#include <glas/bindings/adaptor.hpp>
#include <glas2/concept/coordinate_sparse_matrix.hpp>

namespace glas {

  namespace detail {

    template <typename E>
    typename row_index_array_result_type<E>::type coo_index1_storage( row_orientation, E& e ) {
      return row_index_array( e ) ;
    }

    template <typename E>
    typename column_index_array_result_type<E>::type coo_index1_storage( column_orientation, E& e ) {
      return column_index_array( e ) ;
    }

    template <typename E>
    typename row_index_array_result_type<E>::type coo_index2_storage( column_orientation, E& e ) {
      return row_index_array( e ) ;
    }

    template <typename E>
    typename column_index_array_result_type<E>::type coo_index2_storage( row_orientation, E& e ) {
      return column_index_array( e ) ;
    }

  } // namespace detail


  template <typename E>
  struct CoordinateSparseStructureCollection< E, typename boost::enable_if< boost::mpl::and_< external_type<E>, boost::mpl::not_< boost::is_const<E> > > >::type >
  : boost::mpl::and_< boost::is_same< typename boost::numeric::bindings::detail::property_at< E, boost::numeric::bindings::tag::data_structure >::type, boost::numeric::bindings::tag::coordinate_sparse >
                    , boost::is_same< typename boost::numeric::bindings::detail::property_at< E, boost::numeric::bindings::tag::entity >::type, boost::numeric::bindings::tag::matrix >
                    >
  {} ;


} // namespace glas


namespace boost { namespace numeric { namespace bindings { namespace detail {

  template< typename M, typename Id >
  struct adaptor< M, Id, typename std::enable_if< glas2::is< glas2::CoordinateSparseMatrix, M >::value > >::type >
  {
    typedef typename boost::mpl::if_< glas::is_mutable< Id >
                                    , typename M::value_type
                                    , typename M::value_type const
                                    >::type value_type;
    typedef typename boost::mpl::if_< glas::is_mutable< Id >
                                    , typename M::size_type
                                    , typename M::size_type const
                                    >::type index_type;

    typedef mpl::map<
        mpl::pair< tag::value_type, value_type >,
        mpl::pair< tag::index_type, index_type >,
        mpl::pair< tag::entity, tag::matrix >,
        mpl::pair< tag::size_type<1>, std::ptrdiff_t >,
        mpl::pair< tag::size_type<2>, std::ptrdiff_t >,
        mpl::pair< tag::data_structure, tag::coordinate_sparse >,
        mpl::pair< tag::matrix_type, tag::general,  // At the moment only general nonsymmetric matrices are supported.
        mpl::pair< tag::data_order, tag::row_major, // Assumption although this is not guaranteed
        mpl::pair< tag::index_base, boost::mpl::int_< M::index_base > >
    > property_map ;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.num_rows();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.num_columns();
    }

    static value_type* begin_value( Id& id ) {
        return id.data().ptr() ;
    }

    static value_type* end_value( Id& id ) {
        return id.data().ptr() + id.num_nz() ;
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return glas::stride_1( id ) ;
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return glas::stride_2( id ) ;
    }

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
