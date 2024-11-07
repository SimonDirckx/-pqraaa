//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompany_glasing file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_any_glas_hpp
#define cork_coefficient_matrices_any_glas_hpp

#include <cork/concept/has_type_member.hpp>
#include <cork/utility/ref.hpp>
#include <glas2/matrix.hpp>
#include <glas2/vector.hpp>
#include <glas2/sparse.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  template <typename SequenceMatrices>
  class any_glas
  {
    public:
      typedef typename CORK::deref_type< SequenceMatrices>::type                  sequence_type ;
      typedef typename CORK::deref_type<typename sequence_type::value_type>::type matrix_type ;
      typedef typename matrix_type::value_type                                    value_type ;
      typedef typename matrix_type::size_type                                     size_type ;
      typedef typename sequence_type::size_type                                   grade_type ;

      template <typename T>
      using value_type_for = typename std::common_type< value_type, T >::type ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< value_type, T > > ;

    public:
      any_glas( SequenceMatrices C )
      : C_( C )
      {
        assert(CORK::deref(C_).size()!=0) ;
        assert(CORK::deref(CORK::deref(C_)[0]).num_rows()==CORK::deref(CORK::deref(C_)[0]).num_columns()) ;
      }

    public:
      size_type num_rows() const { return CORK::deref(CORK::deref(C_)[0]).num_rows() ; }
      size_type num_columns() const { return CORK::deref(CORK::deref(C_)[0]).num_columns() ; }

      grade_type num_matrices() const { return CORK::deref(C_).size() ; }

      sequence_type const& sequence() const { return CORK::deref(C_) ; }
      sequence_type& sequence() { return CORK::deref(C_) ; }

    public:
      template <typename M>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, M>::value >::type fill( int i, M m ) const {
        assert( m.num_rows()==num_rows() ) ;
        assert( m.num_columns()==num_columns() ) ;
        m = CORK::deref(CORK::deref(C_)[i]) ;
      }

      auto const& operator()( int i ) const {
        assert( i<num_matrices() ) ;
        return CORK::deref(CORK::deref(C_)[i]) ;
      }

    public:
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        assert( i>=0 && i<num_matrices() ) ;
        w += glas2::multiply( CORK::deref(CORK::deref(C_)[i]), x ) ;
      }

     /**
      * Accumulate: return A + a_0 * C_0 + a_1 * C_1 + .. + a_n * C_n
      * 
      * @param coefs    The seaquence of coefficients a_0..a_n
      * @param A        The Matrix on which the result is accumulated.
      * 
      * @returns        The result matrix 
      */
      template <typename Coefs, typename Matrix>
      typename std::enable_if< glas2::is<glas2::CoordinateSparseMatrix, Matrix>::value>::type accumulate( Coefs const& coefs, Matrix& A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
        assert( coefs.size()==num_matrices() ) ;

        for ( int i=0; i<coefs.size(); ++i ) {
           push_back( A, coefs(i) * CORK::deref(CORK::deref(C_)[i]) ) ;
           A.sort_and_compress() ;
        }
      } // accumulate()

      template <typename Coefs, typename Matrix>
      typename std::enable_if< glas2::is<glas2::DenseMatrix, Matrix>::value>::type accumulate( Coefs const& coefs, Matrix A ) const {
        static_assert( std::is_convertible<typename Coefs::value_type, typename Matrix::value_type>::value, "" ) ;
        assert( coefs.size()==num_matrices() ) ;

        for ( int i=0; i<coefs.size(); ++i ) {
           A += coefs(i) * CORK::deref(CORK::deref(C_)[i]) ;
        }
      } // accumulate()

    private:
      SequenceMatrices C_ ;
  } ; // any_glas

} } // namespace CORK::coefficient_matrices

#endif
