#include <vector>
#include <thread>
#include <iostream>
#include <complex>
#include <cmath>
#include <glas2/vector.hpp>
#include <cork/matrix.hpp>
#include <cork/domain/ellipse.hpp>
#include <cork/domain/line_segment.hpp>
#include <cork/approximation/sv_aaa.hpp>
#include <cork/approximation/sv_aaa_max.hpp>
#include <cork/basis/functions.hpp>
#include <eigen3/Eigen/Dense>
#include <PQRAAA/linalg/truncatedQR.hh>

#include <cork/domain/line_segment.hpp>
#include <cork/eigenvalue_selector/inside_domain.hpp>
#include <cork/nonlinear_matrix.hpp>
#include <cork/eigen_pairs.hpp>
#include <cork/set_of_functions.hpp>
#include <cork/aaa_cork.hpp>
#include <cork/domain/ellipse.hpp>
#include <cork/domain/line_segment.hpp>
#include <cork/approximation/sv_aaa.hpp>
#include <cork/approximation/sv_aaa_max.hpp>
#include <cork/basis/functions.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
using namespace std;
using namespace Eigen;
using namespace boost;
template< typename Tval >
using Mat = Eigen::Matrix<Tval,Eigen::Dynamic,Eigen::Dynamic>;

template< typename Tval >
using Vec = Eigen::Vector<Tval,Eigen::Dynamic>;

using Tval = double;
using CTval = complex<Tval>;

Tval sample(Tval dummy)
{
  static boost::mt19937 rng;
  static boost::normal_distribution<> nd(0.,1.);
  return nd(rng);
}
int main(int argc,char ** argv){
    
    bool qraaa = strcmp(argv[1],"0");
    cout<<"bool = "<<qraaa<<endl;
    
    Eigen::VectorXd tau(10);
    tau<<1.,2.,3.,4.,5.,6.,7.,8.,9,10;
    int L = tau.size();
    int n = 100;
    int nZ = 1000;
    ArrayXd Z = ArrayXd::LinSpaced(nZ,-10., 10.);
    Mat<CTval> F = Mat<CTval>::Zero(nZ,n*n);
    Mat<Tval> A = Mat<Tval>::Zero(L,n*n).unaryExpr(std::ptr_fun<Tval,Tval>(sample));
    Vec<Tval> E = Mat<Tval>::Identity(n,n).reshaped();
    Vec<Tval> nrmVec = Vec<Tval>::Zero(nZ);
    for(int indz=0;indz<nZ;indz++){
        Tval s=Z(indz);
        Vec<complex<Tval>> v = Vec<complex<Tval>>::Zero(n);
        for(int l = 0;l<L;++l){
            v(l)=exp(complex<Tval>(0.,-s*tau(l)))*pow(10.,.5*l)/A.row(l).norm();
        }
        Vec<CTval> w = complex(0.,-s)*E.transpose()-v.transpose()*A;
        nrmVec(indz)=w.norm();
        F.row(indz).noalias() = (w/(w.norm())).transpose();
    }

    auto options = std::make_tuple( CORK::options::number_of_sample_points( nZ )
                                , CORK::options::debug_level( 0 )
                                , CORK::options::max_degree( 80 )
                                , CORK::options::aaa_stop_criterion_max( 1e-8 )
                                ) ;
    glas2::contiguous_vector<Tval,int> Zglas(Z.data(),nZ);

    if(qraaa){
        Linalg::TruncColPivHouseholderQR<Eigen::Matrix<CTval,Eigen::Dynamic,Eigen::Dynamic>> qr(F.rows(), F.cols());
        qr.setThreshold(1e-6);
        auto ticQ = std::chrono::high_resolution_clock::now();
        qr.compute(F);
        Mat<CTval> Q = Mat<Tval>::Identity(F.rows(),qr.rank());
        auto H = qr.householderQ();
        H.setLength(qr.rank());
        Q.applyOnTheLeft(H);
        Mat<CTval> R = qr.matrixR().topLeftCorner(qr.rank(),F.cols()).template triangularView<Eigen::Upper>();
        auto tocQ = std::chrono::high_resolution_clock::now();
        auto elapsedQ = std::chrono::duration_cast<chrono::nanoseconds>(tocQ-ticQ);
        cout<<"qr done in"<<elapsedQ.count()/1e9<<endl;
        cout<<"rank : "<<qr.rank()<<endl;
        cout<<"err = "<<(F*qr.colsPermutation()-Q*R).norm()/F.norm()<<endl;
        Vec<Tval> GAMM(Q.cols());
        for(int k = 0;k<Q.cols();++k){
            Tval rkk = std::abs(R(k,k));
            GAMM(k) = rkk;
        }
        cout<<"GAMM = "<<GAMM<<endl;
        Q.noalias()=Q*GAMM.asDiagonal();
        std::cout<<"NaNs in Q? "<<Q.hasNaN()<<endl;
        R.noalias()=GAMM.asDiagonal().inverse()*R;
        
        glas2::contiguous_matrix<CTval,int,glas2::column_major> Qglas(Q.data(),Q.rows(),Q.cols()); //shallow copy
        auto tic = std::chrono::high_resolution_clock::now();
        auto repr_q = CORK::approximation::SV_AAA( Zglas, Qglas, CORK::approximation::sv_aaa_max(), options ) ;
        auto toc = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
        Tval timeQAAA = elapsed.count()/1e9;
        cout<<"nodes are : "<<repr_q.nodes()<<endl;
        cout<<"degr: "<<repr_q.nodes().size()<<endl;
        cout<<"computed in "<<timeQAAA<<" s"<<endl;
    }else{
        glas2::contiguous_matrix<CTval,int,glas2::column_major> Fglas(F.data(),F.rows(),F.cols()); //shallow copy
        auto tic = std::chrono::high_resolution_clock::now();
        auto repr_f = CORK::approximation::SV_AAA( Zglas, Fglas, CORK::approximation::sv_aaa_max(), options ) ;
        auto toc = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
        Tval timeFAAA = elapsed.count()/1e9;
        cout<<"nodes are : "<<repr_f.nodes()<<endl;
        cout<<"degr: "<<repr_f.nodes().size()<<endl;
        cout<<"computed in "<<timeFAAA<<" s"<<endl;
    }

    cout<<"done"<<endl;

return 0;
}