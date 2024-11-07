//  (C) Copyright Karl Meerbergen 2021.
//
//  Use, modification and distribution are subject to the CORK Software 
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt)
#ifndef cork_options_shift_selection_hpp
#define cork_options_shift_selection_hpp

#include <tuple>
#include <cassert>

namespace CORK { namespace options {

  class shift_selection_key {
    public:
      inline shift_selection_key()
      {}

      template <typename ...Ts>
      inline shift_selection_key( std::tuple<Ts...> const& options )
      {}

    public:
      inline auto value() const { return *this ; }
  } ; // class shift_selection_key


  class shift_selection_eigenvalue_selector
  : public shift_selection_key
  {
    public:
      inline shift_selection_eigenvalue_selector( int number, int multiplicity )
      : multiplicity_( multiplicity )
      , number_( number )
      {}

    public:
      inline int multiplicity() const { return multiplicity_ ; }
      inline int number() const { return number_ ; }

    public:
      inline auto value() const { return *this ; }

    private:
      int multiplicity_ ;
      int number_ ;
  } ; // class shift_selection_eigenvalue_selector


  template <typename T>
  class shift_selection
  : public shift_selection_key
  {
    public:
      typedef T value_type ;

    public:
      shift_selection( int number, int multiplicity )
      : multiplicity_( multiplicity )
      , number_(number)
      , shifts_(0)
      {
        assert( number>0 ) ;
        assert( multiplicity>0 ) ;
      }

    public:
      inline auto value() const { return *this ; }

    public:
      int multiplicity() const { return multiplicity_ ; }
      int number() const { return number_ ; }

      auto const& shifts() const { return shifts_ ; }

      template <typename Shifts>
      void shifts( Shifts const& shifts ) { shifts_.resize(shifts.size()) ; shifts_ = shifts ; }

    private:
      int                                multiplicity_ ;
      int                                number_ ;
      glas2::shared_vector< value_type > shifts_ ;
  } ; // class shift_selection<T>

  template <typename T>
  auto shift_selection_from_single_shift( T const& shift ) {
    shift_selection<T> selection( 1, 1 ) ;
    glas2::static_vector<T,1> shifts ; shifts(0) = shift ;
    selection.shifts( shifts ) ;
    return selection ;
  }

  template <typename T>
  auto shift_selection_from_vector( CORK::vector<T> const& shifts, int multiplicity=20 ) {
    shift_selection<T> selection( shifts.size(), multiplicity ) ;
    selection.shifts( shifts ) ;
    return selection ;
  }

  template <typename Domain>
  auto shift_selection_from_domain( Domain const& domain, int number, int multiplicity=20 ) {
    shift_selection<typename Domain::value_type> selection( number, multiplicity ) ;
    selection.shifts( domain.discretize_coarse( number ) ) ;
    return selection ;
  }

  inline auto shift_selection_from_eigenvalue_selector( int number, int multiplicity=20 ) {
    return shift_selection_eigenvalue_selector( number, multiplicity ) ;
  }

  //
  // Selection of shift_selector
  //

  // default
  template <typename EigenvalueSelector>
  auto produce_shift_selection_split( shift_selection_key const& key, EigenvalueSelector const& eigenvalue_selector ) {
    shift_selection< typename EigenvalueSelector::value_type > selection( 1, 1 ) ;
    selection.shifts( eigenvalue_selector.shifts( selection.number()) ) ;
    return selection ;
  }

  // shift_from_eigenvalue_selector
  template <typename EigenvalueSelector>
  auto produce_shift_selection_split( shift_selection_eigenvalue_selector const& key, EigenvalueSelector const& eigenvalue_selector ) {
    shift_selection< typename EigenvalueSelector::value_type > selection( key.number(), key.multiplicity() ) ;
    selection.shifts( eigenvalue_selector.shifts( selection.number() ) ) ;
    std::cout << "shifts : " << selection.shifts() << std::endl ;
    return selection ;
  }

  // shift_from_vector or shift_from_single_shift
  template <typename T, typename EigenvalueSelector>
  auto produce_shift_selection_split( shift_selection<T> const& key, EigenvalueSelector const& eigenvalue_selector ) {
    return key ;
  }

  template <typename EigenvalueSelector, typename Options>
  auto produce_shift_selection( EigenvalueSelector const& eigenvalue_selector, Options const& options ) {
    auto shiftselect = CORK::options::value_of<shift_selection_key>(options);
    return produce_shift_selection_split( shiftselect, eigenvalue_selector ) ;
  }

} } // CORK::options

#endif
