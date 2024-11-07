#include <arpack/naupd.hpp>
#include <arpack/neupd.hpp>
#include <arpack/standard.hpp>
#include <arpack/copy_op.hpp>
#include <arpack/random_vector_op.hpp>
#include <glas2/vector.hpp>
#include <glas2/bindings/vector.hpp>
#include <glas2/matrix.hpp>
#include <glas2/bindings/matrix.hpp>

template <typename F>
struct arpack_wrapper {
  arpack_wrapper( F const& f )
  : f_( f )
  {}

  template <typename T>
  void operator() ( int n, T* in,T* out ) const {
    glas2::contiguous_vector<T,int> v_in( in, n ) ;
    glas2::contiguous_vector<T,int> v_out( out, n ) ;
    F() ( v_in, v_out ) ;
  }

  F const& f_ ;
} ;

template <typename F>
arpack_wrapper<F> make_arpack_wrapper( F const& f ) {
  return arpack_wrapper<F>( f ) ;
}

// functor for laplacian
struct laplacian {
  template <typename V>
  void operator() ( V in, V out ) const {
    out = glas2::range(1,in.size()+1) * in ;
    out( glas2::range(0,out.size()-1) ) -= in( glas2::range(1,in.size()) ) ;
  }
} ;

int main() {
  int n = 1000 ;
  int ncv = 100 ;
  int nev = 5 ;

  // Workspace
  glas2::shared_vector<double> work( 3*n + 3*ncv*ncv + 6*ncv ) ;
  glas2::shared_vector<double> workev( 3*ncv ) ;

  // Starting vector
  glas2::shared_vector<double> v0( n ) ;

  // Krylov iteration vectors
  glas2::shared_matrix<double> V( n, ncv ) ;

  // Problem setting: Laplacian, standard eigenvalue problem
  auto problem = ARPACK::make_standard( 0.0, make_arpack_wrapper( laplacian() ) ) ;

  // Set data for Arnoldi iteration:
  // NEV : number of eigenavalues wanted
  // 30 restarts
  // tolerance 1e-10
  // LM: Largest magnitude eigenvalues wanted
  //
  glas2::fill( v0, 1.0/std::sqrt(n+0.0) ) ;
  ARPACK::data<double> data( std::string("LM"), nev, 30, 1.e-10, v0, V, work, 1 ) ;

  // Compute Krylov space
  int info = ARPACK::naupd( problem, ARPACK::random_vector_op(), data ) ;
  std::cout << "INFO = " << info << std::endl ;
  if (info!=0) return 1 ;

  // Get eigenvalues and eigenvectors
  glas2::shared_vector< fortran_bool_t > select(data.ncv_) ;
  glas2::shared_vector< int > select_c(0) ;
  glas2::shared_vector< double > dr( nev+1 ) ;
  glas2::shared_vector< double > di( nev+1 ) ;

  /*
  // If also eigenvectors are wanted: store eigenvectors in Z
  glas2::shared_matrix< double > Z( n, nev ) ;
  info = ARPACK::neupd( true, select_c, select, dr, di, Z, problem, data, workev ) ;
  if (info!=0) return 1 ;
  */

  // Only eigenvalues:
  info = ARPACK::neupd( false, select_c, select, dr, di, dr, problem, data, workev ) ;
  if (info!=0) return 1 ;

  std::cout << "Eigenvalues real part : " << dr << std::endl ;
  std::cout << "Eigenvalues imaginary part : " << di << std::endl ;

  if (norm_2( di(glas2::range(0,nev)) ) > 1.e-10 ) return 1 ;
  if (norm_2( dr(glas2::range(0,nev)) - glas2::range(n-nev+1,n+1) ) > 1.e-10 ) return 2 ;

  return 0 ;
}
