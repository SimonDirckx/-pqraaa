// Zifan Liu zifan.liu@cs.kuleuven.be

#ifndef glas2_bsp_vector_container_vector_hpp
#define glas2_bsp_vector_container_vector_hpp

#include <glas2/bsp/vector/concept/bsp_vector.hpp>
#include <glas2/bsp/vector/type/distributed_vector.hpp>
#include <glas2/vector/type/contiguous_vector.hpp>

namespace glas2 { namespace bsp {

  template <typename Distribution, typename T, typename S=std::ptrdiff_t>
  class vector
  : public distributed_vector< glas2::contiguous_vector<T,S>, Distribution >
  {
    public:
      typedef glas2::contiguous_vector<T,S>                     internal_type ;
      typedef distributed_vector< internal_type, Distribution > base_type ;


    public:
      // Constructor and Destructor
      vector( S n, Distribution const& distribution )
      : base_type( internal_type( new T[n], n ), distribution )
	, owner_(true)
      {}


      // Copy reference
      vector( vector const& that )
      : base_type( that )
	, owner_(false)
      {}

      ~vector() {
        if (owner_) delete [] this->local().ptr() ;
      }
/*
    public:
      vector& operator=( vector const& that ) {
        base_type(*this) = that ;
        return *this ;
      }

      vector& operator=( vector&& that ) = delete ;

      template <typename E>
      vector& operator=( E const& that ) {
        base_type(*this) = that ;
        return *this ;
      }
*/
    private:
            bool owner_ ;
  } ;

} 

template <typename D, typename T, typename S>
struct concept<bsp::vector<D,T,S>> : concept<typename bsp::vector<D,T,S>::base_type >
{};

} // namespace glas2::bsp



#endif
