//  (C) Copyright Roel Van Beeumen, Wim Michiels & Karl Meerbergen 2016.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)

#ifndef cork_basis_rational_newton_hpp
#define cork_basis_rational_newton_hpp

#include <cork/basis/iterator.hpp>
#include <cork/concept/has_type_member.hpp>
#include <glas2/vector.hpp>
#include <cassert>
#include <limits>
#include <type_traits>

namespace CORK { namespace basis {

  template <typename Points>
  class rational_newton
  {
    public:
      typedef typename std::decay<Points>::type  points_type ;
      typedef typename points_type::size_type    grade_type ;
      typedef typename points_type::size_type    size_type ;

    private:
      typedef typename points_type::value_type   inner_value_type ;
      typedef decltype(std::abs(inner_value_type()))   real_type ;
      static real_type constexpr infinity = std::numeric_limits<real_type>::infinity() ;

    public:
      template <typename ShiftValueType>
      using value_type_for = typename std::common_type< inner_value_type, ShiftValueType >::type ;

      template <typename T>
      using has_value_type_for = has_type_member< std::common_type< inner_value_type, T > > ;

    public:
      explicit rational_newton( Points nodes, Points poles, Points scaling )
      : nodes_( nodes )
      , poles_( poles )
      , scaling_( scaling )
      {
        assert( nodes_.size()==num_terms() ) ;
        assert( scaling_.size()==num_terms() ) ;
        assert( poles_.size()==num_terms()-1 ) ;
      } // rational_newton

    public:
      size_type num_terms() const {
        return nodes_.size() ;
      } // grade

      size_type grade() const {
        return nodes_.size()-1 ;
      } // grade

    public:
      typename std::add_lvalue_reference< typename std::add_const<points_type>::type >::type nodes() const {
        return nodes_ ;
      } // nodes

      typename std::add_lvalue_reference< typename std::add_const<points_type>::type >::type poles() const {
        return poles_ ;
      } // poles

      typename std::add_lvalue_reference< typename std::add_const<points_type>::type >::type scaling() const {
        return scaling_ ;
      } // poles

      typename std::add_lvalue_reference<Points>::type nodes() { return nodes_ ; }
      typename std::add_lvalue_reference<Points>::type poles() { return poles_ ; }
      typename std::add_lvalue_reference<Points>::type scaling() { return scaling_ ; }

    public:
      template <typename T, typename FunctionValues>
      void evaluate( T const& arg, FunctionValues values ) const {
        assert( values.size() == grade()+1 ) ;
        values(0) = 1./scaling_(0) ;
 /*       if (poles_(poles_.size()-1)!=std::numeric_limits<real_type>::infinity())
          values(0) *= (poles_(poles_.size()-1)-arg) ;*/
        for (typename std::decay<FunctionValues>::type::size_type i=1; i<values.size(); ++i) {
          if (poles_(i-1)==std::numeric_limits<real_type>::infinity())
            values(i) = (arg - nodes_(i-1)) * values(i-1) / scaling_(i) ;
          else
            values(i) = ( (arg - nodes_(i-1)) / (scaling_(i)*(poles_(i-1)-arg)) ) * values(i-1) ;
        }
      } // evaluate

    public:
/*      struct functor_type {
        functor_type( value_type const& arg, points_type const& nodes, points_type const& poles )
        : arg_( arg )
        , nodes_( nodes )
        , poles_( poles )
        , index_( 0 )
        {}

        value_type operator() ( value_type const& v ) {
          ++index_ ;
          return v * ( (arg_ - nodes_(index_-1)) / (1.0 - arg_/poles_(index_-1)) ) ;
        }

        value_type         arg_ ;
        points_type const& nodes_ ;
        points_type const& poles_ ;
        grade_type         index_ ;
      } ;

      typedef CORK::basis::iterator<value_type, functor_type> iterator ;
      iterator evaluate_iterator( value_type const& arg ) const {
        return iterator( functor_type( arg, nodes_, poles_ ) ) ;
      } // evaluate_iterator
*/

    private:
      Points nodes_ ;
      Points poles_ ;
      Points scaling_ ;
  } ; // class national_newton

} } // namespace CORK::basis

#endif
