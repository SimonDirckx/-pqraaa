//  (C) Copyright Karl Meerbergen 2020.
//
//  Use, modification and distribution are subject to the 
//  CORK Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt)

#ifndef cork_eigs_cork_result_hpp
#define cork_eigs_cork_result_hpp

#include <cork/utility/ref.hpp>
#include <type_traits>
#include <memory>

namespace CORK { namespace eigs {

//  template <typename Quadruple, typename Linearization, typename EigSelector, typename Information>
//  struct cork_result {
//    typedef typename std::decay<Quadruple>::type     quadruple_type ;
//    typedef typename std::decay<Linearization>::type linearization_type ;
//    typedef typename std::decay<EigSelector>::type   eigenvalue_selector_type ;
//    typedef typename std::decay<Information>::type   information_type ;
//
//    template <typename Quadruple1, typename Linearization1, typename EigSelector1, typename Information1>
//    cork_result( Quadruple1 quadruple, Linearization1 linearization, EigSelector1 eigenvalue_selector, Information1 information )
//    : quadruple( quadruple )
//    , linearization( linearization )
//    , eigenvalue_selector( eigenvalue_selector )
//    , information( information )
//    {}
//
//    Quadruple     quadruple ;
//    Linearization linearization ;
//    EigSelector   eigenvalue_selector ;
//    Information   information ;
//  } ;

  template <typename Quadruple, typename Linearization, typename EigSelector, typename Information>
  struct cork_result {

    typedef typename deref_type<Quadruple>::type         quadruple_type;
    typedef typename deref_type<Linearization>::type     linearization_type;
    typedef typename deref_type<EigSelector>::type       eig_selector_type;
    typedef typename deref_type<Information>::type       information_type;

    cork_result( Quadruple quadruple, Linearization linearization, EigSelector eigenvalue_selector, Information information) 
    : quadruple_(quadruple)
    , linearization_(linearization)
    , eigenvalue_selector_(eigenvalue_selector)
    , information_(information)
    { }

    public:
    quadruple_type& quadruple() { return CORK::deref(quadruple_) ; }
    linearization_type& linearization() { return CORK::deref(linearization_) ; }
    eig_selector_type& eigenvalue_selector() { return CORK::deref(eigenvalue_selector_) ; }
    information_type& information() { return CORK::deref(information_) ; }

    quadruple_type const& quadruple() const { return deref(quadruple_) ; }
    linearization_type const& linearization() const { return deref(linearization_) ; }
    eig_selector_type const& eigenvalue_selector() const { return deref(eigenvalue_selector_) ; }
    information_type const& information() const { return deref(information_) ; }

    private:
      Quadruple       quadruple_;
      Linearization   linearization_;
      EigSelector     eigenvalue_selector_;
      Information     information_;
  } ; // struct cork_result

  /*template <typename Quadruple, typename Linearization, typename EigSelector, typename Information>
  auto make_cork_result( Quadruple quadruple, Linearization linearization, EigSelector eigenvalue_selector, Information information ) {
    typedef typename std::conditional< std::is_rvalue_reference< Quadruple >::value
                                     , typename std::decay< Quadruple >::type
                                     , Quadruple >::type  quadruple_type ;

    typedef typename std::conditional< std::is_rvalue_reference< Linearization >::value
                                     , typename std::decay< Linearization >::type
                                     , Linearization >::type  linearization_type ;

    typedef typename std::conditional< std::is_rvalue_reference< EigSelector >::value
                                     , typename std::decay< EigSelector >::type
                                     , EigSelector >::type  eigenvalue_selector_type ;

    typedef typename std::conditional< std::is_rvalue_reference< Information >::value
                                     , typename std::decay< Information >::type
                                     , Information >::type  information_type ;

    return cork_result< quadruple_type, linearization_type, eigenvalue_selector_type, information_type >( quadruple, linearization, eigenvalue_selector, information ) ;
  } // make_cork_result()*/

} } // namespace CORK::eigs


#endif
