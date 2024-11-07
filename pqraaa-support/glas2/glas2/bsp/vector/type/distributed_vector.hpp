//Zifan Liu zifan.liu@gmail.com

#ifndef glas2_bsp_vector_type_distributed_vector_hpp
#define glas2_bsp_vector_type_distributed_vector_hpp

#include <glas2/concept/concept.hpp>
//#include <glas2/bsp/backend/default/vector/assign.hpp>
//#include <glas2/bsp/backend/current_backend.hpp>
#include <glas2/bsp/vector/concept/bsp_vector.hpp>

namespace glas2 { namespace bsp {

  template <typename V, typename D>
  class distributed_vector {
    public:
      typedef typename V::value_type value_type ;
      typedef typename V::size_type size_type;

    public:
      distributed_vector( V const& v, D const& d )
      : v_( v )
      , d_( d )
      {}

    public:
       // Copy reference !!
       distributed_vector( distributed_vector const& that )
       : v_(that.v_)
       , d_( that.d_ )
       {}


/*    public:
      distributed_vector operator=( distributed_vector const& that ) {
        assert( d_==that.d_) ;
        v_ = that.v_ ;
        return *this ;
      }

      template <typename E>
      distributed_vector operator=( E const& that ) {
        assign( glas2::current_backend(), *this, that ) ;
        return *this ;
      }
*/
    public:
      typedef V local_type ;
      typedef D distribution_type ;
      local_type  local() const& { return v_ ; } //Can read and write through this function
      distribution_type const& distribution() const { return d_ ; }

    private:
      V        v_ ;
      D const& d_ ;
  } ;

} } // namespace glas::bsp

namespace glas2 {
template <typename V, typename D>
struct concept< bsp::distributed_vector<V,D> >
	: bsp::BSPVector
{} ;

} // namespace glas

#endif
