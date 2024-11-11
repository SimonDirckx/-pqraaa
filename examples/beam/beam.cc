// 
// Project     : PQRAAA
// File        : beam.cc
// Description : test qraaa on clamped sandwich beam reference problem
// Author      : Simon Dirckx
// Copyright   : KU Leuven Dept. CS 2024
//


#include <vector>
#include <thread>
#include <iostream>
#include <complex>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <PQRAAA/sv-aaa/qraaa.hh>

using namespace std;
using namespace Eigen;
using namespace boost;
template< typename Tval >
using Mat = Eigen::Matrix<Tval,Eigen::Dynamic,Eigen::Dynamic>;

template< typename Tval >
using Vec = Eigen::Vector<Tval,Eigen::Dynamic>;

using Tval = double;
using CTval = complex<Tval>;



int main(int argc,char** argv){
    /*
    read (vectorized, conformally-structured) sparse matrices from .dat files
    */
    char filenum[255];
    if(argc==1){
        strcpy(filenum,"168");
    }else{
        strcpy(filenum,argv[1]);
        if(!strcmp(filenum,"168")==0 & !strcmp(filenum,"840")==0 & !strcmp(filenum,"3360")==0){
            cout<<"reached1"<<endl;
            std::__throw_invalid_argument("ERROR: select 168,840 or 3360");
            return 0;
        }
    }

    bool qraaa = (strcmp(argv[2],"1")==0);

    char fmfile[255] = "include/files/beam/vm_";
    char fefile[255] = "include/files/beam/ve_";
    char fvfile[255] = "include/files/beam/vv_";

    strcat(strcat(fmfile,filenum),".dat");
    strcat(strcat(fefile,filenum),".dat");
    strcat(strcat(fvfile,filenum),".dat");
    
     
    double num; int N=0;
    ifstream fm,fe,fv;
    fm.open(fmfile);    while (fm >> num){N++;} fm.close();
    
    Vec<double> vm(N),ve(N),vv(N);
    fm.open(fmfile);    fe.open(fefile);    fv.open(fvfile);
    
    int i = 0;
    while (fm >> num){vm(i++)=num;} i=0;    fm.close();
    while (fe >> num){ve(i++)=num;} i=0;    fe.close();
    while (fv >> num){vv(i++)=num;}         fv.close(); 

    /*
    model parameters and system definition
    
    */

    double p    = .675;
    double G0   = 3.504e5;
    double Ginf = 3.062e9;
    double tau  = 8.230e-9;

    auto g = [p,G0,Ginf,tau] (Tval lam)
    {   CTval z = pow((CTval(0.,lam*tau)),p);
        return (G0 + Ginf*z)/(1. + z);
    };

    /*
    construct F

    */
    int nZ = 1000;
    Eigen::ArrayXd Z = Eigen::ArrayXd::LinSpaced(nZ,200.,30000.);
    Mat<CTval> F = Mat<CTval>::Zero(nZ,N);
    for(int i=0;i<nZ;++i){
        Vec<CTval> w = (ve+g(Z(i))*vv-Z(i)*Z(i)*vm);
        F.row(i).noalias() = w.transpose();
    }
    Vec<Tval> nrmVec(N);
    for(int i=0;i<N;++i){
        nrmVec(i) = F.col(i).norm();
        F.col(i)/=nrmVec(i);
    }
    
    double tol = 1e-8;
    if(qraaa){
        QRAAA::infoType info;
        QRAAA::AAAopts opts;
        opts.tol = tol;
        opts.max_degree = 20;
        auto repr_f=QRAAA::qr_aaa(F,Z,opts,info);
        QRAAA::summarize(info);

        //validate
        int nZtest              = 2513  ;  
        Vec<Tval> err(nZtest)           ;
        Eigen::ArrayXd Ztest    = Eigen::ArrayXd::LinSpaced(nZtest,200.,30000.);
        
        Vec<CTval> ftest(N);
        Vec<CTval> rtest(N);
        
        for(int i=0;i<nZtest;++i){
            Tval z = Ztest(i);
            ftest = (ve+g(z)*vv-z*z*vm);
            
            QRAAA::eval(repr_f,z,rtest);
            for(int idx = 0;idx<N;++idx){rtest(idx)*=nrmVec(idx);}
            cout<<(ftest-rtest).rows()<<","<<(ftest-rtest).cols()<<endl;
            err(i) = (ftest-rtest).array().abs().maxCoeff()/ftest.array().abs().maxCoeff();
        }
        cout<<"Max. err.: "<<err.maxCoeff()<<endl;
        
    }else{

        QRAAA::infoType info;
        QRAAA::AAAopts opts;
        opts.tol = tol;
        opts.max_degree = 20;
        auto repr_f=QRAAA::sv_aaa(F,Z,opts,info);
        QRAAA::summarize(info);
        
        
        //validate
        int nZtest = 2513;  Eigen::ArrayXd Ztest = Eigen::ArrayXd::LinSpaced(nZtest,200.,30000.);
        Vec<CTval> ftest = Vec<CTval>::Zero(N);
        Vec<Tval> err(nZtest);
        Vec<CTval> rtest(N);
        
        for(int i=0;i<nZtest;++i){
            
            Tval z = Ztest(i);
            ftest = (ve+g(z)*vv-z*z*vm);
            
            QRAAA::eval(repr_f,z,rtest);
            for(int idx = 0;idx<N;++idx){rtest(idx)*=nrmVec(idx);}
            err(i) = (ftest-rtest).array().abs().maxCoeff()/ftest.array().abs().maxCoeff();
        
        }
        cout<<"Max. err.: "<<err.maxCoeff()<<endl;
    }

}