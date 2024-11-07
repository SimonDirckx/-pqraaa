//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_basis4cork_union_of_bases_hpp
#define cork_basis4cork_union_of_bases_hpp

#include <cork/basis/union_of_bases.hpp>
#include <cork/basis4cork/union_of_bases_handle.hpp>
#include <cork/basis4cork/basis4cork.hpp>
#include <glas2/matrix.hpp>
#include <glas2/sparse.hpp>
#include <glas2/vector.hpp>
#include <type_traits>
#include <cassert>
#ifndef NDEBUG
#include <limits>
#endif

namespace CORK { namespace basis4CORK {

  template <typename B1, typename B2>
  class basis4CORK< basis::union_of_bases<B1,B2> >
  {
    public:
      template <typename T>
      using value_type_for = typename basis::union_of_bases<B1,B2>::template value_type_for<T> ;

      typedef typename basis::union_of_bases<B1,B2>::size_type              size_type ;

    private:
      typedef ::CORK::basis4CORK::basis4CORK< typename basis::union_of_bases<B1,B2>::basis_1_type > basis4cork_type_1 ;
      typedef ::CORK::basis4CORK::basis4CORK< typename basis::union_of_bases<B1,B2>::basis_2_type > basis4cork_type_2 ;

    public:
      explicit basis4CORK( basis::union_of_bases<B1,B2> const& basis )
      : b1_( basis.basis_1() )
      , b2_( basis.basis_2() )
      , columns_( b2_.size() )
      {
        if (columns_.size()>0) {
          columns_(0) = 0 ;
          columns_(glas2::range_from_end(1,0)) = glas2::range(std::max(1,b1_.size()),size()) ;
        }
      }

    public:
      // Basis4CORK
      size_type num_terms() const { return b1_.num_terms() + b2_.num_terms() - 1 ; }
      size_type size() const { return std::max<size_type>(b1_.size(), std::max<size_type>( b2_.size(), b1_.size()+b2_.size() -1 ) ) ; }

    public:
      //Basis4CORK
      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_M( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        fill( A, 0.0 ) ;
        if (b1_.size()>1) b1_.fill_M( A( glas2::range(0,b1_.size()-1), glas2::range(0,b1_.size())) ) ;
        auto A_2 = A( glas2::range_from_end(std::max(0,b1_.size()-1),0), columns_) ;
        glas2::matrix<typename AA::value_type> temp( A_2.num_rows(), A_2.num_columns() ) ;
        temp = A_2 ;
        b2_.fill_M( temp ) ;
        A_2 = temp ;
      } // fill_M

      template <typename AA>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, AA >::value >::type fill_N( AA A ) const {
        assert( A.num_rows() == size()-1 ) ;
        assert( A.num_columns() == size() ) ;
        fill( A, 0.0 ) ;
        if (b1_.size()>1) b1_.fill_N( A( glas2::range(0,b1_.size()-1), glas2::range(0,b1_.size())) ) ;
        auto A_2 = A( glas2::range_from_end(std::max(0,b1_.size()-1),0), columns_) ;
        glas2::matrix<typename AA::value_type> temp( A_2.num_rows(), A_2.num_columns() ) ;
        temp = A_2 ;
        b2_.fill_N( temp ) ;
        A_2 = temp ;
      } // fill_N

    public:
      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_M( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assert( Z.num_rows()==U.num_rows()-1 ) ;
        assert( Z.num_rows()==size()-1 ) ;
        b2_.multiply_M( U(columns_, glas2::all()), Z(glas2::range_from_end(std::max(0,b1_.size()-1),0), glas2::all()) ) ;
        if (b1_.size()>1) b1_.multiply_M( U(glas2::range(0,b1_.size()), glas2::all()), Z(glas2::range(0,b1_.size()-1), glas2::all()) ) ;
      } // multiply_M()

      // Basis4CORK
      template <typename UM, typename ZM, typename Backend=glas2::default_backend>
      void multiply_N( UM const& U, ZM Z, Backend const& backend=glas2::default_backend() ) const {
        assert( Z.num_rows()==U.num_rows()-1 ) ;
        assert( Z.num_rows()==size()-1 ) ;
        b2_.multiply_N( U(columns_, glas2::all()), Z(glas2::range_from_end(std::max(0,b1_.size()-1),0), glas2::all()) ) ;
        if (b1_.size()>1) b1_.multiply_N( U(glas2::range(0,b1_.size()), glas2::all()), Z(glas2::range(0,b1_.size()-1), glas2::all()) ) ;
      } // multiply_N()

    public:
      template <typename ValueType>
      using handle_type = union_of_bases_handle< typename basis4cork_type_1::template handle_type<ValueType>
                                               , typename basis4cork_type_2::template handle_type<ValueType>
                                               > ;

      template <typename ValueType>
      handle_type<ValueType> handle() const {
        return handle_type<ValueType>( b1_.template handle<ValueType>(), b2_.template handle<ValueType>() ) ;
      }

    private:
      basis4cork_type_1          b1_ ;
      basis4cork_type_2          b2_ ;
      glas2::shared_vector<int>  columns_ ;
  } ; // union_of_bases

} } // namespace CORK::basis4CORK

#endif
