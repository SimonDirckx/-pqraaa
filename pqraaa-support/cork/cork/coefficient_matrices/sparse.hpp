//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_sparse_hpp
#define cork_coefficient_matrices_sparse_hpp

#include <cork/concept/has_type_member.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/sparse.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  namespace detail {
    template <typename SequenceMatrices>
    struct sparse_default_value {
      typedef typename std::remove_const< typename std::remove_reference<SequenceMatrices>::type>::type sequence_type ;
      typedef typename sequence_type::value_type                                                        matrix_type ;
      typedef typename matrix_type::value_type                                                          type ;
    } ;
  }

  template <typename SequenceMatrices>
  class sparse
  {
    public:
      typedef typename std::decay< SequenceMatrices>::type                sequence_type ;
      typedef typename sequence_type::value_type                          matrix_type ;
      typedef typename matrix_type::value_type                            value_type ;
      typedef typename matrix_type::size_type                             size_type ;
      typedef typename sequence_type::size_type                           grade_type ;

      template <typename T>
      using value_type_for = typename std::common_type< value_type, T >::type ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< value_type, T > > ;

    public:
      sparse( SequenceMatrices C )
      : C_( C )
      {
        assert(C_.size()!=0) ;
        assert(C_[0].num_rows()==C_[0].num_columns()) ;
      }

    public:
      size_type num_rows() const { return C_[0].num_rows() ; }
      size_type num_columns() const { return C_[0].num_columns() ; }

      grade_type num_matrices() const { return C_.size() ; }

      sequence_type const& sequence() const { return C_ ; }
      sequence_type& sequence() { return C_ ; }

    public:
      template <typename M>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, M>::value >::type fill( int i, M m ) const {
        assert( m.num_rows()==num_rows() ) ;
        assert( m.num_columns()==num_columns() ) ;
        m = C_[i] ;
      }

      auto const& operator()( int i ) const {
        assert( i<num_matrices() ) ;
        return C_[i] ;
      }

    public:
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        assert( i>=0 && i<num_matrices() ) ;
        w += glas2::multiply( C_[i], x ) ;
      }

      // accumulate
      template <typename Coefs, typename Matrix>
      typename std::enable_if< glas2::is<glas2::CoordinateSparseMatrix, Matrix>::value>::type accumulate( Coefs const& coefs, Matrix& A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
        assert( coefs.size()==num_matrices() ) ;

        for ( int i=0; i<coefs.size(); ++i ) {
           push_back( A, coefs(i) * C_[i] ) ;
           A.sort_and_compress() ;
        }
      } // accumulate()

      template <typename Coefs, typename Matrix>
      typename std::enable_if< glas2::is<glas2::DenseMatrix, Matrix>::value>::type accumulate( Coefs const& coefs, Matrix A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;

        glas2::matrix< typename Matrix::value_type > temp( num_rows(), num_columns() ) ;
        assert( coefs.size()==num_matrices() ) ;

        for ( int i=0; i<coefs.size(); ++i ) {
          temp = C_[i] ;
           A += coefs(i) * temp ;
        }
      } // accumulate()

    private:
      SequenceMatrices C_ ;
  } ; // sparse

} } // namespace CORK::coefficient_matrices

#endif
