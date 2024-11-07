#include <glas2/matrix/concept/dense_matrix.hpp>
#include <cassert>

namespace glas2 {

  template <typename S, typename O>
  class contiguous_iterator {
    public:
      typedef S size_type ;
      typedef O orientation ;

    public:
      contiguous_iterator( size_type size_1, size_type size_2 )
      : size_1_(size_1)
      , size_2_(size_2)
      , index_(0)
      {}

      size_type size_1() { return size_1_ ;}
      size_type size_2() { return size_2_ ;}

      void index_1_pp() {}
      void index_2_pp() { ++ index_ ; }
      size_type operator*() const {return index_ ;}

      size_type stride_1() const {return 1;}
      size_type stride_2() const {return size_1_;}

    private:
      size_type size_1_ ;
      size_type size_2_ ;
      size_type index_ ;
  } ;

}
