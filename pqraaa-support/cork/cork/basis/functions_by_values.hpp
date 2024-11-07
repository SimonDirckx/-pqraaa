#ifndef cork_basis_functions_by_values_hpp
#define cork_basis_functions_by_values_hpp

#include <type_traits>

namespace CORK { namespace basis {

 /**
  * contains a list of functions specified by sampled values.
  *
  * 
  */
  template <typename SampleSequence, typename ValueTable, typename I=typename std::decay< SampleSequence >::type::size_type>
  class functions_by_values
  {
    public:
      typedef typename std::decay< SampleSequence >::type   sample_sequence_type;
      typedef typename std::decay< ValueTable >::type       value_table_type; 
      typedef I                                             size_type;

    public:
      explicit functions_by_values( SampleSequence samples, ValueTable values )
      : samples_( samples )
      , values_( values )
      {}

    public:
     /**
      * Returns the amount of functions stored as values in the table.
      */
      size_type num_terms() const {
        return values_.num_columns();
      }

    public:
     /**
      * Evaluates all the functions in the point arg and stores the result in the vector functionvalues.
      */
      template <typename ValueType, typename FunctionValues> 
      void evaluate( ValueType const& arg, FunctionValues values ) {
        for (typename std::decay<FunctionValues>::type::size_type i=0; i<values.size(); ++i) {
          values(i) = i;
        }
      } // evaluate()

    private:
      SampleSequence samples_;
      ValueTable values_;

  } ; //class functions_by_values

}} // CORK::basis

#endif // cork_basis_functions_by_values_hpp
