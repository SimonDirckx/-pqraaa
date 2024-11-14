// 
// Project     : PQRAAA
// File        : randomDelays.cc
// Description : test qraaa on random delay DE reference problem
// Author      : Simon Dirckx
// Copyright   : KU Leuven Dept. CS 2024
//

#include <vector>
#include <thread>
#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include <PQRAAA/sv-aaa/qraaa.hh>

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
  static boost::uniform_01<> nd;
  return nd(rng);
}
int main(int argc,char ** argv){
    
    int n           = 20;
    bool qraaa      = false;
    int  n_cores    = 1;
    double tol      = 1e-4;
    if (argc>=2){
        n = stoi(argv[1]);
    }
    if(argc>=3){
        bool qraaa = (strcmp(argv[2],"1")==0);
    }
    if(argc>=4){
        n_cores = std::max(1,std::min(8,stoi(argv[3])));
    }
    if(argc>=5){
        tol = stod(argv[4]);
    }
    
    /*
    construct F
    */
    //////////////////////////////////////////////////////////////
    Vec<Tval> tau = Vec<Tval>::Zero(20);
    for(int i=1;i<=tau.size();++i){
        tau(i-1)=1.*i;
    }
    int L = tau.size();
    int nZ = 1000;
    int N = n*n;
    ArrayXd Z = ArrayXd::LinSpaced(nZ,-10., 10.);
    Mat<CTval> F = Mat<CTval>::Zero(nZ,N);
    Mat<Tval> A = Mat<Tval>::Zero(L,N).unaryExpr(std::ptr_fun<Tval,Tval>(sample));
    Vec<Tval> E = Mat<Tval>::Identity(n,n).reshaped();
    Vec<Tval> nrmVec = Vec<Tval>::Zero(nZ);
    for(int indz=0;indz<nZ;indz++){
        Tval s=Z(indz);
        Vec<complex<Tval>> v = Vec<complex<Tval>>::Zero(L);
        for(int l = 0;l<L;++l){
            Tval nrm = A.row(l).reshaped(n,n).cwiseAbs().rowwise().sum().maxCoeff();
            v(l)=exp(complex<Tval>(0.,-s*tau(l)))*pow(10.,.1*l)/nrm;
        }
        Vec<CTval> w = complex(0.,s)*E.transpose()+v.transpose()*A;
        F.row(indz).noalias() = w;
    }
    //////////////////////////////////////////////////////////////


    /*
    approximate & output info
    */
    //////////////////////////////////////////
    QRAAA::infoType info;
    QRAAA::AAAopts  opts;
    opts.tol        = tol;
    opts.max_degree = 100;
    opts.qr         = true;
    opts.n_cores    = n_cores;
    auto repr_f     = QRAAA::sv_aaa(F,Z,opts,info);
    QRAAA::summarize(info);
    //////////////////////////////////////////

    /*
    validate
    */
    /////////////////////////////////////////////////////////////////////////////////////
    int nZtest = 2513;   
    Eigen::ArrayXd Ztest = Eigen::ArrayXd::LinSpaced(nZtest,-10,10);
    Vec<CTval> ftest = Vec<CTval>::Zero(N);
    Vec<CTval> rtest = Vec<CTval>::Zero(N);
    Vec<Tval> err(nZtest);
    for(int j = 0;j<nZtest;++j){
        Tval s=Ztest(j);
        Vec<complex<Tval>> v = Vec<complex<Tval>>::Zero(n);
        for(int l = 0;l<L;++l){
            Tval nrm = A.row(l).reshaped(n,n).cwiseAbs().rowwise().sum().maxCoeff();
            v(l)=exp(complex<Tval>(0.,-s*tau(l)))*pow(10.,.1*l)/nrm;
        }
        Vec<CTval> w = complex(0.,s)*E.transpose()+v.transpose()*A;
        ftest = w;
        QRAAA::eval(repr_f,Ztest(j),rtest);
        err(j) = (ftest-rtest).array().abs().maxCoeff()/ftest.array().abs().maxCoeff();
    }
    cout<<"Maximum error = "<<err.maxCoeff()<<endl;
    /////////////////////////////////////////////////////////////////////////////////////

return 0;
}