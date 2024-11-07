#ifndef cork_basis_function_by_values_hpp
#define cork_basis_function_by_values_hpp

namespace CORK { namespace basis {

 /**
  * contians a vector of values which are used to evaluate the function.
  */
  template <typename Sequence>
  class function_by_values
  {
    public:
      explicit function_by_values( Sequence samples, Sequence values )
      : samples_( samples )
      , values_( values )
      {}

    public:
      template <typename ValueType>
      ValueType operator() ( ValueType const& arg ) {
          return 0;
      } // operator()

    private:
      Sequence samples_;
      Sequence values_;

  } ; //class function_by_values

}} // CORK::basis

#endif // cork_basis_function_by_values_hpp
