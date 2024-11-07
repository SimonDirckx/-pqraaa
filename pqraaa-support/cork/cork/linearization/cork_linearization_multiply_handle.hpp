//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_linearization_cork_linearization_multiply_handle_hpp
#define cork_linearization_cork_linearization_multiply_handle_hpp

#include <cork/linearization/info.hpp>
#include <cork/coefficient_matrices/coefficient_matrices4cork.hpp>
#include <cork/utility/timer.hpp>
#include <cork/utility/value_type.hpp>
#include <cork/utility/value_type_for.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <type_traits>
#include <exception>

namespace CORK { namespace linearization {

  template <typename ValueType, typename Basis4CORK, typename MatrixIterator, typename CoefficientMatrices>
  class cork_linearization_multiply_handle
  {
    public:
      typedef typename std::decay<Basis4CORK>::type                     basis_type ;
      typedef typename std::decay<MatrixIterator>::type                 matrix_iterator_type ;
      typedef typename std::decay<CoefficientMatrices>::type            coefficient_matrices_type ;

      typedef typename basis_type::size_type                            size_type ;
      typedef typename std::common_type< typename coefficient_matrices_type::template value_type_for<ValueType>
                                       , typename basis_type::template value_type_for<ValueType>
                                       >::type                          value_type ;

    private:
      typedef decltype(std::abs(value_type()))                                                         real_value_type ;
      typedef coefficient_matrices::coefficient_matrices4CORK< value_type, coefficient_matrices_type > coefficient_matrices_4cork_type ;

    public:
      cork_linearization_multiply_handle( Basis4CORK basis, MatrixIterator matrix_iterator, CoefficientMatrices coefficient_matrices, info& information )
      : basis_( basis )
      , matrix_iterator_( matrix_iterator )
      , coefficient_matrices_( coefficient_matrices )
      , coefficient_matrices_4cork_( coefficient_matrices_ )
      , information_( information )
      {
        information_.linearization = "CORK Linearization" ;
      }

    public:
      // Multiply with alpha*A + beta*B
      //
      // Input vector is in Compact format: (I\otimes Q) U_in
      // The output vector is
      //  (w_out )
      //  ( 0    ) + (I\otimes Q) U_out
      //  ( |    )
      //  ( 0    )
      // The first row of U_out will be zero.
      //
      template <typename Alpha, typename Beta, typename QQ, typename UUIn, typename WOut, typename UUOut>
      void multiply_AB( Alpha const& alpha, Beta const& beta, QQ const& Q, UUIn const& U_in, WOut w_out, UUOut U_out ) {
        static_assert( std::is_convertible< Alpha, value_type >::value, "CORK::linerization::cork_linearization: Alpha cannot be converted to value_type of output vector" ) ;
        static_assert( std::is_convertible< Beta, value_type >::value, "CORK::linerization::cork_linearization: Beta cannot be converted to value_type of output vector" ) ;
        static_assert( std::is_convertible< typename UUIn::value_type, value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (U_in) cannot be converted to value_type of output vector" ) ;
        static_assert( std::is_convertible< typename QQ::value_type, value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (Q) cannot be converted to value_type of output vector" ) ;
        static_assert( std::is_convertible< value_type, typename WOut::value_type >::value, "CORK::linerization::cork_linearization: value_type cannot be converted to value_type of output vector (w_out)" ) ;

        assert( w_out.size()==Q.num_rows() ) ;
        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;
        assert( U_in.num_columns()==Q.num_columns() ) ;
        assert( U_out.num_columns()==Q.num_columns() ) ;

        glas2::matrix<value_type> U2( U_out.num_rows()-1, U_out.num_columns() ) ;

        // Compute lower blocks

        fill( U_out, 0.0 ) ;

        auto Zr = U_out( glas2::range_from_end(1,0), glas2::all() ) ;
        basis_.multiply_M( alpha*U_in, Zr ) ;
        basis_.multiply_N( beta*U_in, U2 ) ;
        Zr += U2 ;

        // Compute first block.
        glas2::fill( w_out, 0.0 ) ;

        coefficient_matrices_4cork_.initialize_schedule( Q, glas2::range(0, coefficient_matrices_.num_matrices() ) ) ;
        if (alpha!=0.0) matrix_iterator_.schedule_a_0( coefficient_matrices_4cork_/*coefficient_matrices_4cork_.accumulator()*/, [&U_in,&alpha] (auto i) { return alpha*U_in(i, glas2::all()) ; } ) ;
        if (beta!=0.0) matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_/*.accumulator()*/, [&U_in,&beta] (auto i) { return beta*U_in(i, glas2::all()) ; } ) ;
        CORK::timer timer ; timer.tic() ;
        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Q, w_out ) ;
        information_.time_of_matvecs += timer.toc() ;
      } // multiply_AB()

      template <typename QQ, typename UUIn, typename W0, typename UUOut>
      void multiply_A( QQ const& Q, UUIn const& U_in, W0 w0, UUOut U_out ) {
        typedef typename W0::value_type multiply_value_type ;
        static_assert( std::is_convertible< typename UUIn::value_type, multiply_value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (U_in) cannot be converted to value_type of output vector" ) ;
        static_assert( std::is_convertible< typename QQ::value_type, multiply_value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (Q) cannot be converted to value_type of output vector" ) ;

        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        fill( U_out, 0.0 ) ;

        auto Zr = U_out( glas2::range_from_end(1,0), glas2::all() ) ;
        basis_.multiply_M( U_in, Zr ) ;

        // Compute first block.
        glas2::fill( w0, 0.0 ) ;

        coefficient_matrices_4cork_.initialize_schedule( Q, glas2::range(0, coefficient_matrices_.num_matrices() ) ) ;
        matrix_iterator_.schedule_a_0( coefficient_matrices_4cork_/*.accumulator()*/, [&U_in] (auto i) { return U_in(i, glas2::all()) ; } ) ;
        CORK::timer timer ; timer.tic() ;
        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Q, w0 ) ;
        information_.time_of_matvecs += timer.toc() ;
      } // multiply_A()

      template <typename QQ, typename UUIn, typename W0, typename UUOut>
      void multiply_B( QQ const& Q, UUIn const& U_in, W0 w0, UUOut U_out ) {
        typedef typename W0::value_type multiply_value_type ;
        static_assert( std::is_convertible< typename UUIn::value_type, multiply_value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (U_in) cannot be converted to value_type of output vector" ) ;
        static_assert( std::is_convertible< typename QQ::value_type, multiply_value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (Q) cannot be converted to value_type of output vector" ) ;

        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        fill( U_out, 0.0 ) ;

        auto Zr = U_out( glas2::range_from_end(1,0), glas2::all() ) ;
        basis_.multiply_N( U_in, Zr ) ;

        // Compute first block.
        glas2::fill( w0, 0.0 ) ;

        coefficient_matrices_4cork_.initialize_schedule( Q, glas2::range(0, coefficient_matrices_.num_matrices() ) ) ;
        matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_/*.accumulator()*/, [&U_in] (auto i) { return U_in(i, glas2::all()) ; } ) ;
        CORK::timer timer ; timer.tic() ;
        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, Q, w0 ) ;
        information_.time_of_matvecs += timer.toc() ;
      } // multiply_B()

      //
      // Multiplication in full format, i.e. with Q = eye(n,n)
      //
      // U_in and U_out must be column_major
      // U_out(i,j) : i corresponds to degree, j corresponds to dof of the matrix
      //
      template <typename UUIn, typename UUOut>
      void multiply_B( UUIn const& U_in, UUOut U_out ) {
        typedef typename UUOut::value_type multiply_value_type ;
        static_assert( std::is_convertible< typename UUIn::value_type, multiply_value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (U_in) cannot be converted to value_type of output vector" ) ;

        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        fill( U_out, 0.0 ) ;

        auto Zr = U_out( glas2::range_from_end(1,0), glas2::all() ) ;
        auto w0 = U_out( 0, glas2::all() ) ;
        basis_.multiply_N( U_in, Zr ) ;

        // Compute first block.
        coefficient_matrices_4cork_.initialize_schedule( U_out, glas2::range(0, coefficient_matrices_.num_matrices() ) ) ;
        matrix_iterator_.schedule_b_0( coefficient_matrices_4cork_/*.accumulator()*/, [&U_in] (auto i) { return U_in(i,glas2::all()) ; } ) ;
        CORK::timer timer ; timer.tic() ;
        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, w0 ) ;
        information_.time_of_matvecs += timer.toc() ;
      } // multiply_B()

      template <typename UUIn, typename UUOut>
      void multiply_A( UUIn const& U_in, UUOut U_out ) {
        typedef typename UUOut::value_type multiply_value_type ;
        static_assert( std::is_convertible< typename UUIn::value_type, multiply_value_type >::value, "CORK::linerization::cork_linearization: value_type of input vector (U_in) cannot be converted to value_type of output vector" ) ;

        assert( U_in.num_rows()==basis_.size() ) ;
        assert( U_out.num_rows()==basis_.size() ) ;

        fill( U_out, 0.0 ) ;

        auto Zr = U_out( glas2::range_from_end(1,0), glas2::all() ) ;
        auto w0 = U_out( 0, glas2::all() ) ;
        basis_.multiply_M( U_in, Zr ) ;

        // Compute first block.
        coefficient_matrices_4cork_.initialize_schedule( U_out, glas2::range(0, coefficient_matrices_.num_matrices() ) ) ;
        matrix_iterator_.schedule_a_0( coefficient_matrices_4cork_/*.accumulator()*/, [&U_in] (auto i) { return U_in(i, glas2::all()) ; } ) ;
        CORK::timer timer ; timer.tic() ;
        coefficient_matrices_4cork_.apply_scheduled( coefficient_matrices_, w0 ) ;
        information_.time_of_matvecs += timer.toc() ;
      } // multiply_A()

    private:
      Basis4CORK                       basis_ ;
      MatrixIterator                   matrix_iterator_ ;
      CoefficientMatrices              coefficient_matrices_ ;
      coefficient_matrices_4cork_type  coefficient_matrices_4cork_ ;
      info&                            information_ ;
  } ; // default_linearization


} } // namespace CORK::linearization

#endif
