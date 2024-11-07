//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_coefficient_matrices_combined_hpp
#define cork_coefficient_matrices_combined_hpp

#include <cork/concept/has_type_member.hpp>
#include <glas2/vector.hpp>
#include <glas2/matrix.hpp>
#include <cassert>
#include <type_traits>

namespace CORK { namespace coefficient_matrices {

  template <typename CoefficientMatrices, typename Combinations>
  class combined
  {
    public:
      typedef typename std::decay< CoefficientMatrices >::type                                                                         coefficient_matrices_type ;
      typedef typename std::decay< Combinations >::type                                                                                combinations_type ;
      typedef typename coefficient_matrices_type::grade_type                                                                           grade_type ;
      typedef typename coefficient_matrices_type::size_type                                                                            size_type ;

      typedef typename std::common_type< typename coefficient_matrices_type::value_type, typename std::decay<combinations_type>::type::value_type>::type value_type ;

      template <typename T>
      using value_type_for = typename std::common_type< typename coefficient_matrices_type::template value_type_for<T>, typename combinations_type::value_type>::type ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< value_type, T > > ;

    public:
      combined( CoefficientMatrices coefficient, Combinations combinations )
      : coefficient_matrices_( coefficient )
      , combinations_( combinations )
      {
        assert( combinations.num_rows() == coefficient_matrices_.num_matrices() ) ;
      }

    public:
      size_type num_rows() const { return coefficient_matrices_.num_rows() ; }
      size_type num_columns() const { return coefficient_matrices_.num_columns() ; }

      grade_type num_matrices() const { return combinations_.num_columns() ; }

    public:
      coefficient_matrices_type const& coefficient_matrices() const { return coefficient_matrices_ ; }
      combinations_type const& combinations() const { return combinations_ ; }

    private:
      template <typename X, typename EnableIf=void>
      struct multiply_add_traits {
        static auto const& apply( X const& x ) { return x ; }
      } ;

      template <typename V, typename S1, typename S2>
      struct multiply_add_traits< glas2::binary_operation<S1,glas2::binary_operation<S2,V,glas2::multiplies>,glas2::multiplies>, typename std::enable_if< glas2::is<glas2::ContiguousDenseVector,V>::value >::type > {
        template <typename X>
        //static auto apply( X const& x ) { return x;}
        static auto apply( X const& x ) { return (x.scalar()*x.vector().scalar()) * x.vector().vector() ;}
      } ;

    public:
      // Not efficient for a range of i's.
      template <typename X, typename W>
      void multiply_add( grade_type i, X const& x, W w ) const {
        assert( i>=0 && i<num_matrices() ) ;
        glas2::vector<typename W::value_type> xx( x.size() ) ;
        for (size_type j=0; j<combinations_.num_rows(); ++j) {
          //coefficient_matrices_.multiply_add( j, combinations_(j,i)*x, w ) ;
          //xx = copy(combinations_(j,i)*x) ;
          //coefficient_matrices_.multiply_add( j, xx, w ) ;
          auto xc = combinations_(j,i)*x ;
          coefficient_matrices_.multiply_add( j, multiply_add_traits<decltype(xc)>::apply(xc), w ) ;
        }
      }

      template <typename Coefs, typename Matrix>
      void accumulate( Coefs const& coefs, Matrix& matrix ) const {
        coefficient_matrices_.accumulate( multiply(combinations_, coefs), matrix ) ;
      } // accumulate()

    public:
      template <typename M>
      typename std::enable_if< glas2::is< glas2::DenseMatrix, M >::value >::type fill( int i, M m ) const {
        assert( m.num_rows()==num_rows() ) ;
        assert( m.num_columns()==num_columns() ) ;
        assert(i<num_matrices()) ;
        //glas2::matrix<value_type> m_temp(m.num_rows(),m.num_columns()) ;
        glas2::fill( m, 0.0 ) ;
        for (size_type j=0; j<combinations_.num_rows(); ++j) {
          if (combinations_(j,i)!=0.0) {
            //coefficient_matrices_.fill(j, m_temp) ;
            //m += combinations_(j,i) * m_temp ;
            m += combinations_(j,i) * coefficient_matrices_(j) ;
          }
        }
      } // fill()

    private:
      template <typename Selection, typename Lambda, typename W>
      void accumulate( Selection const& selection, Lambda const& lambda, W w ) const {
        assert( w.size()==coefficient_matrices_.num_rows() ) ;
        glas2::vector< typename W::value_type > temp( coefficient_matrices_.num_columns() ) ;

        // Make the sum for each matrix in coefficient_matrices_
        // Lambda my_lambda computes the contribution for each of those matrices.
        coefficient_matrices_.accumulate( glas2::range( 0, combinations_.num_rows() ), [&](auto i) {return this->my_lambda<decltype(i),decltype(temp),Selection,Lambda>(i,temp,selection,lambda);}, w ) ;
      } // accumulate()

    private:
      template <typename I, typename Temp, typename Selection, typename Lambda>
      auto const& my_lambda( I i, Temp& temp, Selection const& selection, Lambda const& lambda ) const {
        fill( temp, 0.0 ) ;
        for (typename Selection::size_type j=0; j<selection.size(); ++j) {
          temp += combinations_(i,selection(j)) * lambda(selection(j)) ;
        }
        return temp ;
      } // lambda

    private:
      CoefficientMatrices coefficient_matrices_ ;
      Combinations        combinations_ ;
  } ; // combined


  template <typename CoefficientMatrices, typename Combinations>
  auto make_combined( CoefficientMatrices const& coefs, Combinations const& combinations ) {
    return combined< CoefficientMatrices, Combinations >( coefs, combinations ) ;
  }


  template <typename CoefficientMatrices>
  struct is_combined
  : std::false_type
  {} ;


  template <typename CoefficientMatrices, typename Combinations>
  struct is_combined< combined< CoefficientMatrices, Combinations > >
  : std::true_type
  {} ;

} } // namespace CORK::coefficient_matrices

#endif
