#ifndef glas2_gnuplot_glas2std_hpp
#define glas2_gnuplot_glas2std_hpp

namespace glas2 { namespace gnuplot {

  namespace detail {

    template <typename X>
    class glas2std_type {
      public:
        glas2std_type( X const& x )
        : x_( x )
        {}

      public:
        typedef typename X::size_type  size_type ;
        typedef typename X::value_type value_type ;

        size_type size() const { return x_.size() ; }
        value_type operator[]( size_type i ) const { return x_(i) ; }

      private:
        X const& x_ ;
    } ;

  } // namespace detail

  template <typename X>
  detail::glas2std_type< X > glas2std( X const& x ) {
    return detail::glas2std_type< X >( x ) ;
  }

} } // namespace glas2::gnuplot


#endif
